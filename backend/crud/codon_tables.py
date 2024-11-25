from uuid import UUID

import sqlalchemy as sa

from ..models import CodonTable
from .base import BaseRepository


class CodonTableRepository(BaseRepository):
    def list_defaults(self, organism: str | None):
        stmt = (
            sa.select(CodonTable)
            .where(CodonTable.user_id == sa.null())
            .order_by(CodonTable.organism, CodonTable.name)
        )

        if organism:
            stmt = stmt.where(CodonTable.organism == organism)

        return self.session.execute(stmt).scalars().all()

    def list_defaults_and_from_user(self, user_id: UUID, organism: str | None):
        """List both user and default tables."""
        stmt = (
            sa.select(CodonTable)
            .where((CodonTable.user_id == user_id) | (CodonTable.user_id == sa.null()))
            .order_by(CodonTable.organism, CodonTable.name)
        )

        if organism:
            stmt = stmt.where(CodonTable.organism == organism)

        return self.session.execute(stmt).scalars().all()

    def count_all(self):
        stmt = sa.select(sa.func.count()).select_from(CodonTable)

        return self.session.execute(stmt).scalar_one()

    def count_defaults(self):
        stmt = (
            sa.select(sa.func.count())
            .select_from(CodonTable)
            .where(CodonTable.user_id is None)
        )

        return self.session.execute(stmt).scalar_one()

    def count_from_user(self, user_id: UUID):
        stmt = (
            sa.select(sa.func.count())
            .select_from(CodonTable)
            .where(CodonTable.user_id == user_id)
        )

        return self.session.execute(stmt).scalar_one()

    def get(self, user_id: UUID, id: UUID):
        stmt = sa.select(CodonTable).where(
            (CodonTable.user_id == user_id) | (CodonTable.user_id == sa.null()),
            CodonTable.id == id,
        )

        return self.session.execute(stmt).scalar_one_or_none()

    def add(self, data: dict) -> UUID | None:
        stmt = sa.insert(CodonTable).values(data).returning(CodonTable.id)
        result = self.execute_with_commit(stmt)

        return result.scalar_one_or_none()

    def update(self, id: UUID, data: dict):
        stmt = sa.update(CodonTable).where(CodonTable.id == id).values(data)
        self.execute_with_commit(stmt)

    def delete(self, user_id: UUID, id: UUID):
        stmt = sa.delete(CodonTable).where(
            CodonTable.user_id == user_id, CodonTable.id == id
        )
        self.execute_with_commit(stmt)
