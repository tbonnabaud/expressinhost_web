from uuid import UUID

import sqlalchemy as sa

from ..models import User
from .base import BaseRepository


class UserRepository(BaseRepository):
    def list(self) -> sa.Sequence[User]:
        stmt = sa.select(User)

        return self.session.execute(stmt).scalars().all()

    def get(self, id: UUID) -> User | None:
        stmt = sa.select(User).where(User.id == id)

        return self.session.execute(stmt).scalar_one_or_none()

    def add(self, data: dict) -> UUID | None:
        stmt = sa.insert(User).values(data).returning(User.id)
        result = self.session.execute(stmt)
        self.session.commit()

        return result.scalar_one_or_none()

    def update(self, id: UUID, data: dict) -> None:
        stmt = sa.update(User).where(User.id == id).values(data)
        self.session.execute(stmt)
        self.session.commit()

    def delete(self, id: UUID) -> None:
        stmt = sa.delete(User).where(User.id == id)
        self.session.execute(stmt)
        self.session.commit()
