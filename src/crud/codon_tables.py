from uuid import UUID

import sqlalchemy as sa

from ..models import CodonTable
from .base import BaseRepository


class CodonTableRepository(BaseRepository):
    def list(self, user_id: UUID | None, organism: str | None):
        stmt = sa.select(CodonTable)

        if user_id is None:
            stmt = stmt.where(CodonTable.user_id == sa.null())

        else:
            stmt = stmt.where(
                (CodonTable.user_id == user_id) | (CodonTable.user_id == sa.null())
            )

        if organism:
            stmt = stmt.where(CodonTable.organism == organism)

        return self.session.execute(stmt).scalars().all()

    def get(self, name: str):
        stmt = sa.select(CodonTable).where(CodonTable.name == name)

        return self.session.execute(stmt).scalar_one_or_none()

    def add(self, data: dict):
        stmt = sa.insert(CodonTable).values(data).returning(CodonTable.name)
        result = self.session.execute(stmt)
        self.session.commit()

        return result.scalar_one_or_none()

    def update(self, name: str, data: dict):
        stmt = sa.update(CodonTable).where(CodonTable.name == name).values(data)
        self.session.execute(stmt)
        self.session.commit()

    def delete(self, name: str):
        stmt = sa.delete(CodonTable).where(CodonTable.name == name)
        self.session.execute(stmt)
        self.session.commit()
