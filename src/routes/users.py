from uuid import UUID

from fastapi import APIRouter, status
from fastapi.exceptions import HTTPException

from ..authentication import hash_password
from ..crud.users import UserRepository
from ..custom_types import SessionDependency
from ..schemas import User, UserForm

router = APIRouter(prefix="/api", tags=["users"])


@router.get("/users", response_model=list[User])
def list_users(session: SessionDependency):
    return UserRepository(session).list()


@router.get("/users/{user_id}", response_model=User)
def get_user(session: SessionDependency, user_id: UUID):
    user = UserRepository(session).get(user_id)

    if user is None:
        raise HTTPException(status.HTTP_404_NOT_FOUND, "Not found.")

    return user


@router.post("/users", response_model=UUID)
def add_user(session: SessionDependency, data: UserForm):
    data.password = hash_password(data.password)

    return UserRepository(session).add(data.model_dump())


@router.put("/users/{user_id}")
def update_user(session: SessionDependency, user_id: UUID, data: UserForm):
    return UserRepository(session).update(user_id, data.model_dump())


@router.delete("/users/{user_id}")
def delete_user(session: SessionDependency, user_id: UUID):
    return UserRepository(session).delete(user_id)
