from typing import Annotated
from uuid import UUID

from fastapi import APIRouter, Depends, status
from fastapi.exceptions import HTTPException
from fastapi.security import OAuth2PasswordRequestForm

from ..authentication import (
    ACCESS_TOKEN_EXPIRE_DELTA,
    check_is_admin,
    check_password,
    create_access_token,
    get_current_user,
    hash_password,
    oauth2_scheme,
)
from ..crud.users import UserRepository
from ..custom_types import SessionDependency
from ..schemas import Token, User, UserForm

router = APIRouter(tags=["Users"])


@router.get("/users", response_model=list[User], dependencies=[Depends(check_is_admin)])
def list_users(session: SessionDependency):
    return UserRepository(session).list()


@router.get("/users/me", response_model=User)
def get_me(session: SessionDependency, token: Annotated[str, Depends(oauth2_scheme)]):
    current_user = get_current_user(session, token)

    if current_user is None:
        raise HTTPException(status.HTTP_404_NOT_FOUND, "Not found.")

    return current_user


@router.get(
    "/users/{user_id}", response_model=User, dependencies=[Depends(check_is_admin)]
)
def get_user(session: SessionDependency, user_id: UUID):
    user = UserRepository(session).get(user_id)

    if user is None:
        raise HTTPException(status.HTTP_404_NOT_FOUND, "Not found.")

    return user


@router.post("/token")
def log_in_user(
    session: SessionDependency,
    form_data: Annotated[OAuth2PasswordRequestForm, Depends()],
):
    user = UserRepository(session).get_by_email(form_data.username)

    if not user or not check_password(form_data.password, user.hashed_password):
        raise HTTPException(status_code=400, detail="Incorrect username or password")

    access_token = create_access_token(
        data={"sub": user.email}, expires_delta=ACCESS_TOKEN_EXPIRE_DELTA
    )

    return Token(access_token=access_token, token_type="bearer")


@router.post("/users", response_model=UUID)
def add_user(session: SessionDependency, data: UserForm):
    user = data.model_dump()
    user["hashed_password"] = hash_password(data.password)
    del user["password"]

    return UserRepository(session).add(user)


@router.put("/users/me")
def update_me(
    session: SessionDependency,
    token: Annotated[str, Depends(oauth2_scheme)],
    data: UserForm,
):
    current_user = get_current_user(session, token)

    updated_user = data.model_dump()
    updated_user["hashed_password"] = hash_password(data.password)
    del updated_user["password"]

    return UserRepository(session).update(current_user.id, updated_user)


@router.delete("/users/me")
def delete_me(
    session: SessionDependency, token: Annotated[str, Depends(oauth2_scheme)]
):
    current_user = get_current_user(session, token)

    return UserRepository(session).delete(current_user.id)


@router.put("/users/{user_id}", dependencies=[Depends(check_is_admin)])
def update_user(session: SessionDependency, user_id: UUID, data: UserForm):
    updated_user = data.model_dump()
    updated_user["hashed_password"] = hash_password(data.password)
    del updated_user["password"]

    return UserRepository(session).update(user_id, updated_user)


@router.delete("/users/{user_id}", dependencies=[Depends(check_is_admin)])
def delete_user(session: SessionDependency, user_id: UUID):
    return UserRepository(session).delete(user_id)
