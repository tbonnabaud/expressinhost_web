import sqlalchemy as sa

from ..models import CodonTranslation
from .base import BaseRepository


class CodonTranslationRepository(BaseRepository):
    def list(self, codon_table_name: str):
        stmt = sa.select(CodonTranslation).where(
            CodonTranslation.codon_table_name == codon_table_name
        )

        return self.session.execute(stmt).scalars().all()

    def get(self, codon_table_name: str, codon: str):
        stmt = sa.select(CodonTranslation).where(
            CodonTranslation.codon_table_name == codon_table_name,
            CodonTranslation.codon == codon,
        )

        return self.session.execute(stmt).scalar_one_or_none()

    def add_batch(self, data_batch: list[dict]):
        stmt = sa.insert(CodonTranslation).values(data_batch)
        self.session.execute(stmt)
        self.session.commit()

    def update(self, codon_table_name: str, codon: str, data: dict):
        stmt = (
            sa.update(CodonTranslation)
            .where(
                CodonTranslation.codon_table_name == codon_table_name,
                CodonTranslation.codon == codon,
            )
            .values(data)
        )
        self.session.execute(stmt)
        self.session.commit()

    def delete(self, codon_table_name: str, codon: str):
        stmt = sa.delete(CodonTranslation).where(
            CodonTranslation.codon_table_name == codon_table_name,
            CodonTranslation.codon == codon,
        )
        self.session.execute(stmt)
        self.session.commit()
