from uuid import UUID

from fastapi import APIRouter

from ..authentication import OptionalTokenDependency, TokenDependency, get_current_user
from ..crud.codon_tables import CodonTableRepository
from ..crud.codon_translations import CodonTranslationRepository
from ..database import SessionDependency
from ..schemas import CodonTable, UserCodonTableFormWithTranslations, CodonTranslation

router = APIRouter(tags=["Codon tables"])


def assign_codon_table_id(codon_table_id: UUID, tr: CodonTranslation) -> dict:
    """
    The function `assign_codon_table_id` assigns a UUID codon table ID to a CodonTranslation object and
    returns the updated dictionary representation of the object.
    """
    tr_dict = tr.model_dump()
    tr_dict["codon_table_id"] = codon_table_id

    return tr_dict


@router.get("/codon-tables", response_model=list[CodonTable])
def list_codon_tables(
    session: SessionDependency,
    token: OptionalTokenDependency,
    organism: str | None = None,
):
    current_user = get_current_user(session, token) if token else None

    if current_user:
        return CodonTableRepository(session).list_defaults_and_from_user(
            current_user.id, organism
        )
    else:
        return CodonTableRepository(session).list_defaults(organism)


@router.get("/users/me/codon-tables/{codon_table_id}", response_model=CodonTable)
def get_user_codon_table(
    session: SessionDependency, token: TokenDependency, codon_table_id: UUID
):
    current_user = get_current_user(session, token)

    return CodonTableRepository(session).get(current_user.id, codon_table_id)


@router.get(
    "/users/me/codon-tables/{codon_table_id}/translations",
    response_model=list[CodonTranslation],
)
def get_user_codon_table_translations(
    session: SessionDependency, token: TokenDependency, codon_table_id: UUID
):
    current_user = get_current_user(session, token)
    codon_table = CodonTableRepository(session).get(current_user.id, codon_table_id)

    if codon_table:
        return CodonTranslationRepository(session).list_from_table(codon_table.id)


@router.post("/users/me/codon-tables")
def add_user_codon_table(
    session: SessionDependency,
    token: TokenDependency,
    data: UserCodonTableFormWithTranslations,
):
    current_user = get_current_user(session, token)

    table = data.model_dump(exclude={"translations"})
    table["user_id"] = current_user.id

    codon_table_id = CodonTableRepository(session).add(table)

    CodonTranslationRepository(session).add_batch(
        list(map(lambda x: assign_codon_table_id(codon_table_id, x), data.translations))
    )

    return codon_table_id


@router.put("/users/me/codon-tables/{codon_table_id}")
def update_user_table_codon_translations(
    session: SessionDependency,
    token: TokenDependency,
    codon_table_id: UUID,
    data: UserCodonTableFormWithTranslations,
):
    current_user = get_current_user(session, token)

    table = data.model_dump(exclude={"translations"})

    codon_table_repo = CodonTableRepository(session)
    codon_table = codon_table_repo.get(current_user.id, codon_table_id)

    if codon_table:
        # Update metadata
        codon_table_repo.update(codon_table.id, table)

        # We replace all old rows by the new ones
        codon_translation_repo = CodonTranslationRepository(session)
        codon_translation_repo.delete_batch(codon_table_id)
        codon_translation_repo.add_batch(
            list(
                map(
                    lambda x: assign_codon_table_id(codon_table_id, x),
                    data.translations,
                )
            )
        )


@router.delete("/users/me/codon-tables/{codon_table_id}")
def delete_user_codon_table(
    session: SessionDependency, token: TokenDependency, codon_table_id: UUID
):
    current_user = get_current_user(session, token)

    return CodonTableRepository(session).delete(current_user.id, codon_table_id)
