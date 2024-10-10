from datetime import UTC, datetime
from uuid import UUID, uuid4

import sqlalchemy as sa
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column


class Base(DeclarativeBase):
    pass


class User(Base):
    __tablename__ = "users"

    id: Mapped[UUID] = mapped_column(
        sa.UUID(as_uuid=True), primary_key=True, default=uuid4
    )
    creation_date: Mapped[datetime] = mapped_column(
        sa.DateTime, default=lambda: datetime.now(UTC)
    )
    email: Mapped[str] = mapped_column(sa.String)
    password: Mapped[str] = mapped_column(sa.String)


class CodonTable(Base):
    __tablename__ = "codon_tables"

    name: Mapped[str] = mapped_column(sa.String, primary_key=True)
    user_id: Mapped[UUID | None] = mapped_column(
        sa.UUID(as_uuid=True),
        sa.ForeignKey("users.id", onupdate="CASCADE", ondelete="CASCADE"),
        nullable=True,
    )
    creation_date: Mapped[datetime] = mapped_column(
        sa.DateTime, default=lambda: datetime.now(UTC)
    )
    organism: Mapped[str] = mapped_column(sa.String)
    custom: Mapped[bool] = mapped_column(sa.Boolean)


class CodonTranslation(Base):
    __tablename__ = "codon_translations"

    codon_table_name: Mapped[str] = mapped_column(sa.String)
    codon: Mapped[str] = mapped_column(sa.String(length=3))
    anti_codon: Mapped[str] = mapped_column(sa.String(length=3))
    amino_acid: Mapped[str] = mapped_column(sa.String(length=3))
    trna_gcn: Mapped[float] = mapped_column(sa.Float)
    corresponding_codon: Mapped[str] = mapped_column(sa.String(length=3))
    wobble_rate: Mapped[float] = mapped_column(sa.Float)

    __table_args__ = (
        sa.PrimaryKeyConstraint(
            "codon_table_name", "codon", name="pk_codon_translations"
        ),
    )


class Result(Base):
    __tablename__ = "results"

    id: Mapped[UUID] = mapped_column(
        sa.UUID(as_uuid=True), primary_key=True, default=uuid4
    )
    user_id: Mapped[UUID | None] = mapped_column(
        sa.UUID(as_uuid=True),
        sa.ForeignKey("users.id", onupdate="CASCADE", ondelete="CASCADE"),
        nullable=True,
    )
    creation_date: Mapped[datetime] = mapped_column(
        sa.DateTime, default=lambda: datetime.now(UTC)
    )
    host_organism: Mapped[str] = mapped_column(sa.String)
    mode: Mapped[str] = mapped_column(sa.String)
    parameters: Mapped[dict] = mapped_column(sa.JSON)


class TunedSequence(Base):
    __tablename__ = "tuned_sequences"

    id: Mapped[UUID] = mapped_column(
        sa.UUID(as_uuid=True), primary_key=True, default=uuid4
    )
    result_id: Mapped[UUID] = mapped_column(
        sa.UUID(as_uuid=True),
        sa.ForeignKey("results.id", onupdate="CASCADE", ondelete="CASCADE"),
    )
    name: Mapped[str] = mapped_column(sa.String)
    input: Mapped[str] = mapped_column(sa.Text)
    output: Mapped[str] = mapped_column(sa.Text)
    percentage_identity: Mapped[float] = mapped_column(sa.Float)
