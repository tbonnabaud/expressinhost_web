from fastapi import APIRouter

from ..authentication import TokenDependency, OptionalTokenDependency, get_current_user
from ..crud.codon_tables import CodonTableRepository
from ..database import SessionDependency
from ..schemas import CodonTable, CodonTableForm

router = APIRouter(tags=["Codon tables"])


@router.get("/codon-tables", response_model=list[CodonTable])
def list_codon_tables(session: SessionDependency, token: OptionalTokenDependency, organism: str | None = None):
    current_user = get_current_user(session, token) if token else None

    if current_user:
        return CodonTableRepository(session).list_defaults_and_from_user(
        current_user.id, organism
    )
    else:
        return CodonTableRepository(session).list_defaults(organism)


# @router.get("/users/me/codon-tables", response_model=list[CodonTable])
# def list_default_and_user_codon_tables(
#     session: SessionDependency,
#     token: TokenDependency,
#     organism: str | None = None,
# ):
#     current_user = get_current_user(session, token)

#     return CodonTableRepository(session).list_defaults_and_from_user(
#         current_user.id, organism
#     )


@router.get("/users/me/codon-tables/{codon_table_name}", response_model=CodonTable)
def get_codon_table(
    session: SessionDependency, token: TokenDependency, codon_table_name: str
):
    current_user = get_current_user(session, token)

    return CodonTableRepository(session).get(current_user.id, codon_table_name)


@router.post("/users/me/codon-tables")
def add_codon_table(
    session: SessionDependency, token: TokenDependency, data: CodonTableForm
):
    current_user = get_current_user(session, token)

    table = data.model_dump()
    table["user_id"] = current_user.id

    return CodonTableRepository(session).add(data)


@router.delete("/users/me/codon-tables/{codon_table_name}")
def delete_codon_table(
    session: SessionDependency, token: TokenDependency, codon_table_name: str
):
    current_user = get_current_user(session, token)

    return CodonTableRepository(session).delete(current_user.id, codon_table_name)
