from datetime import timedelta
from typing import Annotated

from fastapi import APIRouter, Depends, Request, Response
from fastapi.exceptions import HTTPException
from fastapi.responses import RedirectResponse
from fastapi.security import OAuth2PasswordRequestForm

from ..authentication import (
    check_password,
    create_access_token,
    create_refresh_token,
    revoke_refresh_token,
    verify_refresh_token,
)
from ..crud.users import UserRepository
from ..database import SessionDependency
from ..email_service import send_email
from ..settings import settings

router = APIRouter(prefix="/auth", tags=["Authentication"])

ACCESS_TOKEN_EXPIRE_DELTA = timedelta(minutes=15)
REFRESH_TOKEN_EXPIRE_DELTA = timedelta(days=7)
RESET_TOKEN_EXPIRE_DELTA = timedelta(minutes=15)


@router.post("/token")
def get_access_token(
    session: SessionDependency,
    form_data: Annotated[OAuth2PasswordRequestForm, Depends()],
):
    user = UserRepository(session).get_by_email(form_data.username.lower())

    if not user or not check_password(form_data.password, user.hashed_password):
        raise HTTPException(status_code=400, detail="Incorrect username or password")

    access_token = create_access_token(
        user_id=user.id,
        user_role=user.role,
        expires_delta=ACCESS_TOKEN_EXPIRE_DELTA,
    )

    # Create refresh token and store in Redis
    refresh_token = create_refresh_token(
        user_id=user.id, expires_delta=REFRESH_TOKEN_EXPIRE_DELTA
    )

    # Create redirect response to home page
    response = RedirectResponse(url="/", status_code=302)

    # Set httpOnly cookie with JWT access token
    response.set_cookie(
        key="access_token",
        value=access_token,
        httponly=True,
        secure=settings.COOKIE_SECURE,
        samesite="lax",
        max_age=int(ACCESS_TOKEN_EXPIRE_DELTA.total_seconds()),
        path="/",
    )

    # Set httpOnly cookie with refresh token
    response.set_cookie(
        key="refresh_token",
        value=refresh_token,
        httponly=True,
        secure=settings.COOKIE_SECURE,
        samesite="lax",
        max_age=int(REFRESH_TOKEN_EXPIRE_DELTA.total_seconds()),
        path="/api/auth",
    )

    return response


@router.post("/refresh")
def refresh_tokens(request: Request, session: SessionDependency):
    """
    Refresh endpoint that generates new access and refresh tokens.

    Args:
        request: FastAPI Request object to extract refresh token from cookie
        session: Database session

    Returns:
        Response with new access and refresh tokens as httpOnly cookies

    Raises:
        HTTPException: 401 if refresh token is invalid or expired
    """
    # Get refresh token from cookie
    refresh_token = request.cookies.get("refresh_token")

    if not refresh_token:
        raise HTTPException(
            status_code=401,
            detail="Refresh token missing",
        )

    # Verify refresh token and get user ID
    user_id = verify_refresh_token(refresh_token)

    if not user_id:
        raise HTTPException(
            status_code=401,
            detail="Invalid or expired refresh token",
        )

    # Get user from database
    user = UserRepository(session).get(user_id)

    if not user:
        raise HTTPException(
            status_code=401,
            detail="User not found",
        )

    # Revoke old refresh token
    revoke_refresh_token(refresh_token)

    # Create new access token
    new_access_token = create_access_token(
        user_id=user.id,
        user_role=user.role,
        expires_delta=ACCESS_TOKEN_EXPIRE_DELTA,
    )

    # Create new refresh token
    new_refresh_token = create_refresh_token(
        user_id=user.id, expires_delta=REFRESH_TOKEN_EXPIRE_DELTA
    )

    # Create response
    response = Response(status_code=200)

    # Set new access token cookie
    response.set_cookie(
        key="access_token",
        value=new_access_token,
        httponly=True,
        secure=settings.COOKIE_SECURE,
        samesite="lax",
        max_age=int(ACCESS_TOKEN_EXPIRE_DELTA.total_seconds()),
        path="/",
    )

    # Set new refresh token cookie
    response.set_cookie(
        key="refresh_token",
        value=new_refresh_token,
        httponly=True,
        secure=settings.COOKIE_SECURE,
        samesite="lax",
        max_age=int(REFRESH_TOKEN_EXPIRE_DELTA.total_seconds()),
        path="/api/auth",
    )

    return response


@router.post("/logout")
def logout(request: Request):
    """
    Logout endpoint that clears the authentication cookies and revokes refresh token.

    Args:
        request: FastAPI Request object to extract refresh token from cookie

    Returns:
        Response with deleted cookies
    """
    # Get refresh token from cookie and revoke it
    refresh_token = request.cookies.get("refresh_token")
    if refresh_token:
        revoke_refresh_token(refresh_token)

    # Create response and delete both cookies
    response = Response(status_code=200)
    response.delete_cookie(
        key="access_token",
        path="/",
        samesite="lax",
    )
    response.delete_cookie(
        key="refresh_token",
        path="/api/auth",
        samesite="lax",
    )
    return response


@router.get("/password-forgotten")
def send_reset_password_email(
    request: Request, session: SessionDependency, user_email: str
):
    user = UserRepository(session).get_by_email(user_email.lower())

    if user:
        reset_token = create_access_token(
            user_id=user.id,
            user_role=user.role,
            expires_delta=RESET_TOKEN_EXPIRE_DELTA,
            for_reset=True,
        )

        reset_url = f"{request.base_url}reset-password/{reset_token}"

        send_email(
            user.email,
            "Reset password link for ExpressInHost",
            f"Reset password link for ExpressInHost:\n{reset_url}",
        )


@router.get("/reset-password")
def get_reset_data(reset_token: str):
    return {"reset_token": reset_token}
