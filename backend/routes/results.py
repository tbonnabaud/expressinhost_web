from uuid import UUID

from fastapi import APIRouter

from ..authentication import TokenDependency, get_current_user
from ..crud.results import ResultRepository
from ..database import SessionDependency
from ..schemas import Result
from .common import FilterParamDependency

router = APIRouter(tags=["Results"])


@router.get("/users/me/results", response_model=list[Result])
def list_user_results(
    session: SessionDependency,
    token: TokenDependency,
    filter_params: FilterParamDependency,
):
    current_user = get_current_user(session, token)

    return ResultRepository(session).list_from_user(
        current_user.id, filter_params.offset, filter_params.limit
    )


@router.get("/users/me/results/count", response_model=int)
def count_user_results(session: SessionDependency, token: TokenDependency):
    current_user = get_current_user(session, token)

    return ResultRepository(session).count_from_user(current_user.id)


@router.get("/users/me/results/{result_id}", response_model=Result)
def get_user_result(
    session: SessionDependency, token: TokenDependency, result_id: UUID
):
    current_user = get_current_user(session, token)
    return ResultRepository(session).get(current_user.id, result_id)


@router.delete("/users/me/results/{result_id}")
def delete_user_result(
    session: SessionDependency, token: TokenDependency, result_id: UUID
):
    current_user = get_current_user(session, token)
    return ResultRepository(session).delete(current_user.id, result_id)
