from uuid import UUID

from fastapi import APIRouter

from ..authentication import JWTDependency
from ..crud.results import ResultRepository
from ..database import SessionDependency, SessionWithCommitDependency
from ..schemas import Result
from .common import FilterParamDependency

router = APIRouter(tags=["Results"])


@router.get("/users/me/results", response_model=list[Result])
def list_user_results(
    session: SessionDependency,
    jwt: JWTDependency,
    filter_params: FilterParamDependency,
):
    return ResultRepository(session).list_from_user(
        jwt.sub, filter_params.offset, filter_params.limit
    )


@router.get("/users/me/results/count", response_model=int)
def count_user_results(session: SessionDependency, jwt: JWTDependency):
    return ResultRepository(session).count_from_user(jwt.sub)


@router.get("/users/me/results/{result_id}", response_model=Result)
def get_user_result(session: SessionDependency, jwt: JWTDependency, result_id: UUID):
    return ResultRepository(session).get(jwt.sub, result_id)


@router.delete("/users/me/results/{result_id}")
def delete_user_result(
    session: SessionWithCommitDependency, jwt: JWTDependency, result_id: UUID
):
    return ResultRepository(session).delete(jwt.sub, result_id)


@router.delete("/users/me/results")
def delete_user_results(session: SessionWithCommitDependency, jwt: JWTDependency):
    return ResultRepository(session).delete_all(jwt.sub)
