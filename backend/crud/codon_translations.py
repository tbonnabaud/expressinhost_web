from uuid import UUID

import sqlalchemy as sa

from ..models import CodonTranslation
from .base import BaseRepository


class CodonTranslationRepository(BaseRepository):
    def list_from_table(self, codon_table_id: UUID):
        stmt = sa.select(CodonTranslation).where(
            CodonTranslation.codon_table_id == codon_table_id
        )

        return self.session.execute(stmt).scalars().all()

    def get(self, codon_table_id: UUID, codon: str):
        stmt = sa.select(CodonTranslation).where(
            CodonTranslation.codon_table_id == codon_table_id,
            CodonTranslation.codon == codon,
        )

        return self.session.execute(stmt).scalar_one_or_none()

    def add_batch(self, data_batch: list[dict]):
        stmt = sa.insert(CodonTranslation).values(data_batch)
        self.execute_with_commit(stmt)

    def update(self, codon_table_id: UUID, codon: str, data: dict):
        stmt = (
            sa.update(CodonTranslation)
            .where(
                CodonTranslation.codon_table_id == codon_table_id,
                CodonTranslation.codon == codon,
            )
            .values(data)
        )
        self.execute_with_commit(stmt)

    def delete_batch(self, codon_table_id: UUID):
        stmt = sa.delete(CodonTranslation).where(
            CodonTranslation.codon_table_id == codon_table_id
        )
        self.execute_with_commit(stmt)
