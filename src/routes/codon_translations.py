from fastapi import APIRouter

from ..crud.codon_translations import CodonTranslationRepository
from ..custom_types import SessionDependency
from ..schemas import CodonTranslation

router = APIRouter(prefix="/api", tags=["Codon translations"])


@router.get(
    "/codon-tables/{codon_table_name}/codon-translations",
    response_model=list[CodonTranslation],
)
def list_table_codon_translations(session: SessionDependency, codon_table_name: str):
    return CodonTranslationRepository(session).list_from_table(codon_table_name)


@router.get(
    "/codon-tables/{codon_table_name}/codon-translations/{codon}",
    response_model=CodonTranslation,
)
def get_codon_translation(
    session: SessionDependency, codon_table_name: str, codon: str
):
    return CodonTranslationRepository(session).get(codon_table_name, codon)


@router.post("/codon-tables/{codon_table_name}/codon-translations")
def add_table_codon_translations(
    session: SessionDependency,
    data_batch: list[CodonTranslation],
):
    return CodonTranslationRepository(session).add_batch(data_batch)


@router.put("/codon-tables/{codon_table_name}/codon-translations")
def update_table_codon_translations(
    session: SessionDependency,
    codon_table_name: str,
    data_batch: list[CodonTranslation],
):
    repo = CodonTranslationRepository(session)
    # We replace all old rows by the new ones
    repo.delete_batch(codon_table_name)
    repo.add_batch(data_batch)
