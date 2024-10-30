from uuid import UUID

import sqlalchemy as sa

from ..models import User
from .base import BaseRepository


class UserRepository(BaseRepository):
    def list(self):
        stmt = sa.select(User)

        return self.session.execute(stmt).scalars().all()

    def get(self, id: UUID):
        stmt = sa.select(User).where(User.id == id)

        return self.session.execute(stmt).scalar_one_or_none()

    def get_by_email(self, email: str):
        stmt = sa.select(User).where(User.email == email)

        return self.session.execute(stmt).scalar_one_or_none()

    def add(self, data: dict) -> UUID | None:
        stmt = sa.insert(User).values(data).returning(User.id)
        result = self.execute_with_commit(stmt)

        return result.scalar_one_or_none()

    def update(self, id: UUID, data: dict):
        stmt = sa.update(User).where(User.id == id).values(data)
        self.execute_with_commit(stmt)

    def delete(self, id: UUID):
        stmt = sa.delete(User).where(User.id == id)
        self.execute_with_commit(stmt)
