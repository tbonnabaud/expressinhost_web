from uuid import UUID

import sqlalchemy as sa

from ..models import TunedSequence
from .base import BaseRepository


class TunedSequenceRepository(BaseRepository):
    def list(self, result_id: UUID) -> sa.Sequence[TunedSequence]:
        stmt = sa.select(TunedSequence).where(TunedSequence.result_id == result_id)

        return self.session.execute(stmt).scalars().all()

    def get(self, id: UUID) -> TunedSequence | None:
        stmt = sa.select(TunedSequence).where(TunedSequence.id == id)

        return self.session.execute(stmt).scalar_one_or_none()

    def add_batch(self, data_batch: list[dict]) -> None:
        stmt = sa.insert(TunedSequence).values(data_batch)
        self.session.execute(stmt)
        self.session.commit()

    def update(self, id: UUID, data: dict) -> None:
        stmt = sa.update(TunedSequence).where(TunedSequence.id == id).values(data)
        self.session.execute(stmt)
        self.session.commit()

    def delete(self, id: UUID) -> None:
        stmt = sa.delete(TunedSequence).where(TunedSequence.id == id)
        self.session.execute(stmt)
        self.session.commit()
