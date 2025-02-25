from datetime import date, timedelta

from fastapi import APIRouter, Depends

from ..authentication import check_is_admin
from ..crud.run_infos import RunInfoRepository
from ..database import SessionDependency
from ..schemas import RunInfo
from .common import FilterParamDependency

router = APIRouter(tags=["Run infos"], dependencies=[Depends(check_is_admin)])


@router.get("/run-infos", response_model=list[RunInfo])
def get_run_infos(
    session: SessionDependency,
    filter_params: FilterParamDependency,
):
    return RunInfoRepository(session).list(filter_params.offset, filter_params.limit)


@router.get("/run-infos/duration-statistics", response_model=dict[str, timedelta])
def compute_duration_statistics(session: SessionDependency):
    return RunInfoRepository(session).compute_duration_statistics()


@router.get("/run-infos/mode-distribution", response_model=dict[str, int])
def compute_mode_distribution(session: SessionDependency):
    return RunInfoRepository(session).compute_mode_distribution()


@router.get("/run-infos/count-per-day", response_model=dict[date, int])
def compute_run_count_per_day(session: SessionDependency):
    return RunInfoRepository(session).compute_run_count_per_day()


@router.get("/run-infos/sequence-number-statistics", response_model=dict[str, int])
def compute_sequence_number_statistics(session: SessionDependency):
    return RunInfoRepository(session).compute_sequence_number_statistics()
