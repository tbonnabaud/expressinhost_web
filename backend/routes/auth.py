from datetime import timedelta
from typing import Annotated

from fastapi import APIRouter, Depends, Request
from fastapi.exceptions import HTTPException
from fastapi.security import OAuth2PasswordRequestForm

from ..authentication import check_password, create_token
from ..crud.users import UserRepository
from ..database import SessionDependency
from ..schemas import Token

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

    return Token(access_token=access_token, token_type="bearer")


@router.post("/password-forgotten")
def send_reset_password_email(
    request: Request, session: SessionDependency, user_email: str
):
    user = UserRepository(session).get_by_email(user_email.lower())

    if user:
        reset_token = create_token(
            data={"sub": user.email, "for_reset": True},
            expires_delta=RESET_TOKEN_EXPIRE_DELTA,
        )

        return {"reset_token": reset_token, "base_url": str(request.base_url)}

    return {}


@router.get("/reset-password")
def get_reset_data(reset_token: str):
    return {"reset_token": reset_token}
