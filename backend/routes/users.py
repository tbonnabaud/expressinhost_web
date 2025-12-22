from uuid import UUID

from fastapi import APIRouter, Depends, status
from fastapi.exceptions import HTTPException

from ..authentication import (
    TokenDependency,
    check_is_admin,
    get_current_user,
    hash_password,
)
from ..crud.users import UserRepository
from ..database import SessionDependency, SessionWithCommitDependency
from ..schemas import User, UserForm, UserPasswordForm, UserProfileForm, UserRoleForm
from .common import FilterParamDependency

router = APIRouter(tags=["Users"])


@router.get("/users", response_model=list[User], dependencies=[Depends(check_is_admin)])
def list_users(session: SessionDependency, filter_params: FilterParamDependency):
    return UserRepository(session).list(filter_params.offset, filter_params.limit)


@router.get("/users/me", response_model=User)
def get_me(session: SessionDependency, token: TokenDependency):
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


@router.post("/users", response_model=UUID)
def add_user(session: SessionWithCommitDependency, data: UserForm):
    user = data.model_dump()
    user["hashed_password"] = hash_password(data.password)
    del user["password"]

    return UserRepository(session).add(user)


@router.put("/users/me/profile")
def update_me_profile(
    session: SessionWithCommitDependency,
    token: TokenDependency,
    data: UserProfileForm,
):
    current_user = get_current_user(session, token)
    updated_user = data.model_dump()

    return UserRepository(session).update(current_user.id, updated_user)


@router.put("/users/me/password")
def update_me_password(
    session: SessionWithCommitDependency, data: UserPasswordForm, reset: bool = False
):
    current_user = get_current_user(session, data.reset_token, check_for_reset=reset)
    updated_user = {}
    updated_user["hashed_password"] = hash_password(data.password)

    return UserRepository(session).update(current_user.id, updated_user)


@router.delete("/users/me")
def delete_me(session: SessionWithCommitDependency, token: TokenDependency):
    current_user = get_current_user(session, token)

    return UserRepository(session).delete(current_user.id)


@router.put("/users/{user_id}/role", dependencies=[Depends(check_is_admin)])
def update_user_role(
    session: SessionWithCommitDependency, user_id: UUID, data: UserRoleForm
):
    updated_user = data.model_dump()

    return UserRepository(session).update(user_id, updated_user)


@router.delete("/users/{user_id}", dependencies=[Depends(check_is_admin)])
def delete_user(session: SessionWithCommitDependency, user_id: UUID):
    return UserRepository(session).delete(user_id)
