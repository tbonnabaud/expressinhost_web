from datetime import datetime, timedelta, timezone
from typing import Annotated

import bcrypt
import jwt
from fastapi import Depends, status
from fastapi.exceptions import HTTPException
from fastapi.security import OAuth2PasswordBearer

from .crud.users import User, UserRepository
from .database import Session, SessionDependency
from .settings import settings

oauth2_scheme = OAuth2PasswordBearer(tokenUrl="auth/token")
optional_oauth2_scheme = OAuth2PasswordBearer(tokenUrl="auth/token", auto_error=False)

TokenDependency = Annotated[str, Depends(oauth2_scheme)]
OptionalTokenDependency = Annotated[str | None, Depends(optional_oauth2_scheme)]


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


def create_token(data: dict, expires_delta: timedelta | None = None) -> str:
    """
    The `create_token` function generates a JWT token with optional expiration time.

    Args:
        data (dict): The `data` parameter is a dictionary containing the information that you want to
        encode into the token. This information could include user details, permissions, or any other
        data that needs to be included in the token.

        expires_delta (timedelta | None): The `expires_delta` parameter is used to specify the duration
        for which the token will be valid. If a value is provided for `expires_delta`, the token
        will expire after that duration. Default is 15 minutes.

    Returns:
        str: The function `create_token` returns a string which is the encoded JSON Web Token (JWT)
        containing the data provided along with an expiration time.
    """
    to_encode = data.copy()

    if expires_delta:
        expire = datetime.now(timezone.utc) + expires_delta

    else:
        expire = datetime.now(timezone.utc) + timedelta(minutes=15)

    to_encode.update({"exp": expire})
    encoded_jwt = jwt.encode(
        to_encode, settings.JWT_SECRET_KEY, algorithm=settings.JWT_ALGORITHM
    )

    return encoded_jwt


def decode_access_token(token: str) -> dict:
    """
    The function `decode_access_token` decodes a JWT access token using a secret key and algorithm
    specified in the settings.

    Args:
        token (str): JWT token.

    Raises:
        HTTPException: Error 401 Unauthorized.

    Returns:
        dict: A dictionary containing the decoded information from the access token
    """
    try:
        return jwt.decode(
            token, settings.JWT_SECRET_KEY, algorithms=[settings.JWT_ALGORITHM]
        )

    except jwt.InvalidTokenError:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid authentication credentials",
            headers={"WWW-Authenticate": "Bearer"},
        )


def get_current_user(
    session: Session, token: str, check_for_reset: bool = False
) -> User:
    """
    The function `get_current_user` retrieves the current user based on the provided access token and
    session, handling authentication errors if necessary.
    """
    payload = decode_access_token(token)
    email: str | None = payload.get("sub")
    for_reset: bool | None = payload.get("for_reset")

    if email and (check_for_reset and for_reset or not check_for_reset):
        user = UserRepository(session).get_by_email(email)

        if user:
            return user

    raise HTTPException(
        status_code=status.HTTP_401_UNAUTHORIZED,
        detail="Invalid authentication credentials",
        headers={"WWW-Authenticate": "Bearer"},
    )


def check_is_member(session: SessionDependency, token: TokenDependency) -> bool:
    user = get_current_user(session, token)

    if user is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found",
            headers={"WWW-Authenticate": "Bearer"},
        )

    return True


def check_is_admin(session: SessionDependency, token: TokenDependency) -> bool:
    user = get_current_user(session, token)

    if user is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found",
            headers={"WWW-Authenticate": "Bearer"},
        )

    elif user.role != "admin":
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Not admin",
            headers={"WWW-Authenticate": "Bearer"},
        )

    return True
