from uuid import UUID

from fastapi import APIRouter

from ..crud.codon_tables import CodonTableRepository
from ..custom_types import SessionDependency
from ..schemas import CodonTable, UserCodonTableForm

router = APIRouter(prefix="/api", tags=["Codon tables"])


@router.get("/codon-tables", response_model=list[CodonTable])
def list_default_codon_tables(session: SessionDependency, organism: str | None = None):
    return CodonTableRepository(session).list_defaults(organism)


@router.get("/users/{user_id}/codon-tables", response_model=list[CodonTable])
def list_user_codon_tables(
    session: SessionDependency,
    user_id: UUID,
    organism: str | None = None,
):
    return CodonTableRepository(session).list_defaults_and_from_user(user_id, organism)


@router.get("/codon-tables/{codon_table_name}", response_model=CodonTable)
def get_codon_table(session: SessionDependency, codon_table_name: str):
    return CodonTableRepository(session).get(codon_table_name)


@router.post("/codon-tables")
def add_codon_table(session: SessionDependency, data: UserCodonTableForm):
    return CodonTableRepository(session).add(data)


@router.delete("/codon-tables/{codon_table_name}")
def delete_codon_table(session: SessionDependency, codon_table_name: str):
    return CodonTableRepository(session).delete(codon_table_name)
