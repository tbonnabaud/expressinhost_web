from uuid import UUID

import sqlalchemy as sa

from ..models import RunInfo
from .base import BaseRepository


class RunInfoRepository(BaseRepository):
    def list(self, offset: int | None, limit: int | None):
        stmt = (
            sa.select(RunInfo)
            .offset(offset)
            .limit(limit)
            .order_by(sa.desc(RunInfo.creation_date))
        )

        return self.session.execute(stmt).scalars().all()

    def count(self):
        stmt = sa.select(sa.func.count()).select_from(RunInfo)

        return self.session.execute(stmt).scalar_one()

    def get(self, id: UUID):
        stmt = sa.select(RunInfo).where(RunInfo.id == id)

        return self.session.execute(stmt).scalar_one_or_none()

    def add(self, data: dict) -> UUID | None:
        stmt = sa.insert(RunInfo).values(data).returning(RunInfo.id)
        result = self.session.execute(stmt)

        return result.scalar_one_or_none()
