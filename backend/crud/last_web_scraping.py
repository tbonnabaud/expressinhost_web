from datetime import UTC, datetime

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import insert as pg_insert

from ..models import LastWebScraping
from .base import BaseRepository


class LastWebScrapingRepository(BaseRepository):
    def list(self):
        stmt = sa.select(LastWebScraping)

        return self.session.execute(stmt).scalars().all()

    def get_last_from(self, source: str):
        stmt = (
            sa.select(LastWebScraping)
            .where(LastWebScraping.source == source)
            .order_by(LastWebScraping.scraping_date.desc())
        )

        return self.session.execute(stmt).scalar_one_or_none()

    # def add(self, data: dict):
    #     stmt = sa.insert(LastWebScraping).values(data)
    #     self.session.execute(stmt)

    # def delete(self, source: str, release: str):
    #     stmt = sa.delete(LastWebScraping).where(
    #         LastWebScraping.source == source, LastWebScraping.release == release
    #     )
    #     self.session.execute(stmt)

    def upsert(self, source: str, release: str):
        stmt = (
            pg_insert(LastWebScraping)
            .values(source=source, release=release, scraping_date=datetime.now(UTC))
            .on_conflict_do_update(
                index_elements=["source", "release"],
                set_=dict(source=source, release=release),
            )
        )

        self.session.execute(stmt)
