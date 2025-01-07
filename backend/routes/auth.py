from typing import Annotated

from fastapi import APIRouter, Depends
from fastapi.exceptions import HTTPException
from fastapi.security import OAuth2PasswordRequestForm

from ..authentication import (
    ACCESS_TOKEN_EXPIRE_DELTA,
    check_password,
    create_access_token,
)
from ..crud.users import UserRepository
from ..database import SessionDependency
from ..schemas import Token

router = APIRouter(prefix="/auth", tags=["Authentication"])


@router.post("/token")
def get_token(
    session: SessionDependency,
    form_data: Annotated[OAuth2PasswordRequestForm, Depends()],
):
    user = UserRepository(session).get_by_email(form_data.username.lower())

    if not user or not check_password(form_data.password, user.hashed_password):
        raise HTTPException(status_code=400, detail="Incorrect username or password")

    access_token = create_access_token(
        data={"sub": user.email}, expires_delta=ACCESS_TOKEN_EXPIRE_DELTA
    )

    return Token(access_token=access_token, token_type="bearer")
