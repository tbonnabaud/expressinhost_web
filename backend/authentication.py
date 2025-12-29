from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from typing import Annotated
from uuid import UUID, uuid4

import bcrypt
import jwt
from fastapi import Cookie, Depends, status
from fastapi.exceptions import HTTPException

from .redis_client import redis_client
from .settings import settings


@dataclass
class JWTPayload:
    """
    Dataclass representing the JWT token payload structure.

    Attributes:
        sub: User ID (subject)
        role: User's role (e.g., 'admin', 'member')
        exp: Token expiration timestamp
        for_reset: Optional flag indicating if token is for password reset
    """

    sub: UUID
    role: str
    exp: datetime
    for_reset: bool = False


def hash_password(password: str) -> str:
    """
    The function `hash_password` takes a password as input, hashes it using bcrypt, and returns the
    hashed password as a string.
    """
    return bcrypt.hashpw(password.encode(), bcrypt.gensalt()).decode()


def check_password(password: str, hashed_password: str) -> bool:
    """
    The function `check_password` compares a plain text password with a hashed password using bcrypt to
    determine if they match.
    """
    return bcrypt.checkpw(password.encode(), hashed_password.encode())


# Access token
def create_access_token(
    user_id: UUID,
    user_role: str,
    expires_delta: timedelta,
    for_reset: bool = False,
) -> str:
    """
    The `create_access_token` function generates a JWT token with the specified expiration time.

    Args:
        user_id (UUID): The user's unique identifier (stored as 'sub' in the token).
        user_role (str): The user's role (e.g., 'admin', 'member').
        expires_delta (timedelta): The duration for which the token will be valid.
        for_reset (bool): Whether this token is for password reset. Defaults to False.

    Returns:
        str: The encoded JSON Web Token (JWT) containing the user information and expiration time.
    """
    payload = {
        "sub": str(user_id),
        "role": user_role,
    }

    if for_reset:
        payload["for_reset"] = True

    expire = datetime.now(timezone.utc) + expires_delta

    payload["exp"] = expire
    encoded_jwt = jwt.encode(
        payload, settings.JWT_SECRET_KEY, algorithm=settings.JWT_ALGORITHM
    )

    return encoded_jwt


def decode_access_token(token: str) -> JWTPayload:
    """
    The function `decode_access_token` decodes a JWT access token using a secret key and algorithm
    specified in the settings.

    Args:
        token (str): JWT token.

    Raises:
        HTTPException: Error 401 Unauthorized.

    Returns:
        JWTPayload: A dataclass containing the decoded information from the access token
    """
    try:
        payload: dict = jwt.decode(
            token, settings.JWT_SECRET_KEY, algorithms=[settings.JWT_ALGORITHM]
        )
        return JWTPayload(
            sub=UUID(payload["sub"]),
            role=payload["role"],
            exp=datetime.fromtimestamp(payload["exp"], tz=timezone.utc),
            for_reset=payload.get("for_reset", False),
        )

    except jwt.ExpiredSignatureError:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Token has expired",
        )

    except jwt.InvalidTokenError:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid authentication credentials",
        )


# Refresh token
def create_refresh_token(user_id: UUID, expires_delta: timedelta) -> str:
    """
    Create a refresh token (UUID) and store it in Redis with the user's ID.

    Args:
        user_id: The user's ID
        expires_delta: Time delta for token expiration

    Returns:
        str: The generated refresh token (UUID)
    """
    refresh_token = str(uuid4())
    # Store in Redis with user_id as value and expiration time
    redis_client.setex(
        f"refresh_token:{refresh_token}",
        int(expires_delta.total_seconds()),
        str(user_id),
    )
    return refresh_token


def verify_refresh_token(refresh_token: str) -> UUID | None:
    """
    Verify a refresh token and return the associated user ID.

    Args:
        refresh_token: The refresh token to verify

    Returns:
        UUID | None: The user's ID if token is valid, None otherwise
    """
    user_id = redis_client.get(f"refresh_token:{refresh_token}")

    if user_id:
        return UUID(user_id)

    return None


def revoke_refresh_token(refresh_token: str) -> None:
    """
    Revoke a refresh token by deleting it from Redis.

    Args:
        refresh_token: The refresh token to revoke
    """
    redis_client.delete(f"refresh_token:{refresh_token}")


# Dependencies
def extract_jwt_payload(
    access_token: Annotated[str | None, Cookie()] = None,
) -> JWTPayload:
    print("coucou")
    if access_token is None:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Access token missing",
        )

    return decode_access_token(access_token)


def extract_optional_jwt_payload(
    access_token: Annotated[str | None, Cookie()] = None,
) -> JWTPayload | None:
    if access_token is None:
        return None

    return decode_access_token(access_token)


# Aliases
JWTDependency = Annotated[JWTPayload, Depends(extract_jwt_payload)]
OptionalJWTDependency = Annotated[
    JWTPayload | None, Depends(extract_optional_jwt_payload)
]


def check_is_admin(jwt: JWTDependency) -> bool:
    if jwt.role != "admin":
        raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Not admin")

    return True
