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

    def get_by_email(self, email: UUID):
        stmt = sa.select(User).where(User.email == email)

        return self.session.execute(stmt).scalar_one_or_none()

    def add(self, data: dict):
        with self.session.begin():
            stmt = sa.insert(User).values(data).returning(User.id)
            result = self.session.execute(stmt)

            return result.scalar_one_or_none()

    def update(self, id: UUID, data: dict):
        with self.session.begin():
            stmt = sa.update(User).where(User.id == id).values(data)
            self.session.execute(stmt)

    def delete(self, id: UUID):
        with self.session.begin():
            stmt = sa.delete(User).where(User.id == id)
            self.session.execute(stmt)
