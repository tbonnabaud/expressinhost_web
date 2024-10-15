from uuid import UUID

from fastapi import APIRouter

from ..crud.tuned_sequences import TunedSequenceRepository
from ..custom_types import SessionDependency
from ..schemas import TunedSequence

router = APIRouter(tags=["Tuned sequences"])


@router.get("/results/{result_id}/tuned-sequences", response_model=list[TunedSequence])
def list_result_tuned_sequences(session: SessionDependency, result_id: UUID):
    return TunedSequenceRepository(session).list_from_result(result_id)
