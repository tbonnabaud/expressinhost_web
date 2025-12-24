from datetime import timedelta
from typing import Annotated

from fastapi import APIRouter, Depends, Request, Response
from fastapi.exceptions import HTTPException
from fastapi.responses import RedirectResponse
from fastapi.security import OAuth2PasswordRequestForm

from ..authentication import check_password, create_token
from ..crud.users import UserRepository
from ..database import SessionDependency
from ..email_service import send_email
from ..settings import settings

router = APIRouter(prefix="/auth", tags=["Authentication"])

ACCESS_TOKEN_EXPIRE_DELTA = timedelta(hours=12)
RESET_TOKEN_EXPIRE_DELTA = timedelta(minutes=15)


@router.post("/token")
def get_access_token(
    session: SessionDependency,
    form_data: Annotated[OAuth2PasswordRequestForm, Depends()],
):
    user = UserRepository(session).get_by_email(form_data.username.lower())

    if not user or not check_password(form_data.password, user.hashed_password):
        raise HTTPException(status_code=400, detail="Incorrect username or password")

    access_token = create_token(
        data={"sub": user.email}, expires_delta=ACCESS_TOKEN_EXPIRE_DELTA
    )

    # Create redirect response to home page
    response = RedirectResponse(url="/", status_code=302)

    # Set httpOnly cookie with JWT token
    response.set_cookie(
        key="access_token",
        value=access_token,
        httponly=True,
        secure=settings.COOKIE_SECURE,
        samesite="lax",
        max_age=43200,  # 12 hours in seconds (matches ACCESS_TOKEN_EXPIRE_DELTA)
        path="/",
    )

    return response


@router.post("/logout")
def logout():
    """
    Logout endpoint that clears the authentication cookie.

    Returns:
        Response with deleted cookie
    """
    response = Response(status_code=200)
    response.delete_cookie(
        key="access_token",
        path="/",
        samesite="lax",
    )
    return response


@router.get("/password-forgotten")
def send_reset_password_email(
    request: Request, session: SessionDependency, user_email: str
):
    user = UserRepository(session).get_by_email(user_email.lower())

    if user:
        reset_token = create_token(
            data={"sub": user.email, "for_reset": True},
            expires_delta=RESET_TOKEN_EXPIRE_DELTA,
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
