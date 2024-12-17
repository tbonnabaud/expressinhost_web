from uuid import UUID

import sqlalchemy as sa

from ..models import Result
from .base import BaseRepository


class ResultRepository(BaseRepository):
    def list_from_user(self, user_id: UUID, offset: int | None, limit: int | None):
        stmt = (
            sa.select(Result)
            .where(Result.user_id == user_id)
            .offset(offset)
            .limit(limit)
            .order_by(sa.desc(Result.creation_date))
        )

        return self.session.execute(stmt).scalars().all()

    def count_from_user(self, user_id: UUID):
        stmt = (
            sa.select(sa.func.count())
            .select_from(Result)
            .where(Result.user_id == user_id)
        )

        return self.session.execute(stmt).scalar_one()

    def get(self, user_id: UUID, id: UUID):
        stmt = sa.select(Result).where(Result.user_id == user_id, Result.id == id)

        return self.session.execute(stmt).scalar_one_or_none()

    def add(self, data: dict) -> UUID | None:
        stmt = sa.insert(Result).values(data).returning(Result.id)
        result = self.execute_with_commit(stmt)

        return result.scalar_one_or_none()

    def update(self, id: UUID, data: dict):
        stmt = sa.update(Result).where(Result.id == id).values(data)
        self.execute_with_commit(stmt)

    def delete(self, user_id: UUID, id: UUID):
        stmt = sa.delete(Result).where(Result.user_id == user_id, Result.id == id)
        self.execute_with_commit(stmt)

    def delete_all(self, user_id: UUID):
        stmt = sa.delete(Result).where(Result.user_id == user_id)
        self.execute_with_commit(stmt)
