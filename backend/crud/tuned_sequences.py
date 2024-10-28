from uuid import UUID

import sqlalchemy as sa

from ..models import TunedSequence
from .base import BaseRepository


class TunedSequenceRepository(BaseRepository):
    def list_from_result(self, result_id: UUID):
        stmt = sa.select(TunedSequence).where(TunedSequence.result_id == result_id)

        return self.session.execute(stmt).scalars().all()

    def get(self, id: UUID):
        stmt = sa.select(TunedSequence).where(TunedSequence.id == id)

        return self.session.execute(stmt).scalar_one_or_none()

    def add_batch(self, data_batch: list[dict]):
        stmt = sa.insert(TunedSequence).values(data_batch)
        self.execute_with_commit(stmt)

    def update(self, id: UUID, data: dict):
        stmt = sa.update(TunedSequence).where(TunedSequence.id == id).values(data)
        self.execute_with_commit(stmt)

    def delete(self, id: UUID):
        stmt = sa.delete(TunedSequence).where(TunedSequence.id == id)
        self.execute_with_commit(stmt)
