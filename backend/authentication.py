from typing import Annotated

import bcrypt
from fastapi import Depends, status
from fastapi.exceptions import HTTPException
from fastapi.security import OAuth2PasswordBearer

from .crud.users import UserRepository
from .custom_types import Session, SessionDependency

oauth2_scheme = OAuth2PasswordBearer(tokenUrl="token")


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


def get_current_user(session: Session, token: str):
    user = UserRepository(session).get_by_email(token)

    if not user:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid authentication credentials",
            headers={"WWW-Authenticate": "Bearer"},
        )

    return user


def check_is_member(
    session: SessionDependency, token: Annotated[str, Depends(oauth2_scheme)]
):
    return get_current_user(session, token)


def check_is_admin(
    session: SessionDependency, token: Annotated[str, Depends(oauth2_scheme)]
):
    user = get_current_user(session, token)

    if user.role != "admin":
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Not admin",
            headers={"WWW-Authenticate": "Bearer"},
        )

    return user
