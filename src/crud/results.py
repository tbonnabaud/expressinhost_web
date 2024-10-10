from uuid import UUID

import sqlalchemy as sa

from ..models import Result
from .base import BaseRepository


class ResultRepository(BaseRepository):
    def list(self, user_id: UUID) -> sa.Sequence[Result]:
        stmt = sa.select(Result).where(Result.user_id == user_id)

        return self.session.execute(stmt).scalars().all()

    def get(self, id: UUID) -> Result | None:
        stmt = sa.select(Result).where(Result.id == id)

        return self.session.execute(stmt).scalar_one_or_none()

    def add(self, data: dict) -> UUID | None:
        stmt = sa.insert(Result).values(data).returning(Result.id)
        result = self.session.execute(stmt)
        self.session.commit()

        return result.scalar_one_or_none()

    def update(self, id: UUID, data: dict) -> None:
        stmt = sa.update(Result).where(Result.id == id).values(data)
        self.session.execute(stmt)
        self.session.commit()

    def delete(self, id: UUID) -> None:
        stmt = sa.delete(Result).where(Result.id == id)
        self.session.execute(stmt)
        self.session.commit()
