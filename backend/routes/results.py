from uuid import UUID

from fastapi import APIRouter

from ..crud.results import ResultRepository
from ..database import SessionDependency
from ..schemas import Result
from .common import FilterParamDependency

router = APIRouter(tags=["Results"])


@router.get("/users/{user_id}/results", response_model=list[Result])
def list_user_results(
    session: SessionDependency,
    filter_params: FilterParamDependency,
    user_id: UUID,
):
    return ResultRepository(session).list_from_user(
        user_id, filter_params.offset, filter_params.limit
    )


@router.get("/users/{user_id}/results/count", response_model=list[Result])
def count_user_results(
    session: SessionDependency,
    user_id: UUID,
):
    return ResultRepository(session).count_from_user(user_id)


@router.get("/results/{result_id}", response_model=Result)
def get_result(session: SessionDependency, result_id: UUID):
    return ResultRepository(session).get(result_id)


@router.delete("/results/{result_id}")
def delete_result(session: SessionDependency, result_id: UUID):
    return ResultRepository(session).delete(result_id)
