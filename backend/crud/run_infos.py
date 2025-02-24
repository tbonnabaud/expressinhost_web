from datetime import timedelta
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

    def compute_duration_statistics(self) -> dict[str, timedelta]:
        stmt = sa.select(
            sa.func.min(RunInfo.duration).label("min_duration"),
            sa.func.avg(RunInfo.duration).label("avg_duration"),
            sa.func.max(RunInfo.duration).label("max_duration"),
        )

        result = self.session.execute(stmt).one()

        return {
            "min_duration": result.min_duration,
            "avg_duration": result.avg_duration,
            "max_duration": result.max_duration,
        }

    def compute_mode_distribution(self) -> dict[str, int]:
        stmt = sa.select(
            RunInfo.mode, sa.func.count(RunInfo.mode).label("count")
        ).group_by(RunInfo.mode)

        result = self.session.execute(stmt).fetchall()
        return {row.mode: row.count for row in result}

    def compute_run_count_per_day(self) -> dict[sa.Date, int]:
        stmt = (
            sa.select(
                sa.cast(RunInfo.creation_date, sa.Date).label("date"),
                sa.func.count(RunInfo.id).label("count"),
            )
            .group_by(sa.cast(RunInfo.creation_date, sa.Date))
            .order_by("date")
        )

        result = self.session.execute(stmt).fetchall()
        return {row.date: row.count for row in result}
