from datetime import UTC, datetime
from uuid import UUID, uuid4

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship


class Base(DeclarativeBase):
    pass


class User(Base):
    __tablename__ = "users"

    id: Mapped[UUID] = mapped_column(
        sa.UUID(as_uuid=True), primary_key=True, default=uuid4
    )
    creation_date: Mapped[datetime] = mapped_column(
        sa.DateTime(timezone=True), default=lambda: datetime.now(UTC)
    )
    email: Mapped[str] = mapped_column(sa.String, unique=True, index=True)
    hashed_password: Mapped[str] = mapped_column(sa.String)
    role: Mapped[str] = mapped_column(sa.String, default="member")
    full_name: Mapped[str] = mapped_column(sa.String)
    contact_consent: Mapped[bool] = mapped_column(sa.Boolean, server_default="false")


class CodonTable(Base):
    __tablename__ = "codon_tables"

    id: Mapped[UUID] = mapped_column(
        sa.UUID(as_uuid=True), primary_key=True, default=uuid4
    )
    user_id: Mapped[UUID | None] = mapped_column(
        sa.UUID(as_uuid=True),
        sa.ForeignKey("users.id", onupdate="CASCADE", ondelete="CASCADE"),
        nullable=True,
    )
    creation_date: Mapped[datetime] = mapped_column(
        sa.DateTime(timezone=True), default=lambda: datetime.now(UTC)
    )
    name: Mapped[str] = mapped_column(sa.String)
    organism: Mapped[str] = mapped_column(sa.String)

    __table_args__ = (sa.UniqueConstraint("user_id", "name", "organism"),)


class CodonTranslation(Base):
    __tablename__ = "codon_translations"

    codon_table_id: Mapped[UUID] = mapped_column(
        sa.UUID(as_uuid=True),
        sa.ForeignKey("codon_tables.id", onupdate="CASCADE", ondelete="CASCADE"),
        primary_key=True,
    )
    codon: Mapped[str] = mapped_column(sa.String(length=3), primary_key=True)
    anticodon: Mapped[str] = mapped_column(sa.String(length=3))
    amino_acid: Mapped[str] = mapped_column(sa.String(length=3))
    trna_gcn: Mapped[float] = mapped_column(sa.Float)
    wobble_codon: Mapped[str] = mapped_column(sa.String(length=3))
    wobble_rate: Mapped[float] = mapped_column(sa.Float)


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
        sa.DateTime(timezone=True), default=lambda: datetime.now(UTC)
    )
    host_codon_table_id: Mapped[UUID] = mapped_column(
        sa.UUID(as_uuid=True),
        sa.ForeignKey("codon_tables.id", onupdate="CASCADE", ondelete="CASCADE"),
    )
    sequences_native_codon_tables: Mapped[dict] = mapped_column(JSONB)
    mode: Mapped[str] = mapped_column(sa.String)
    slow_speed_threshold: Mapped[float] = mapped_column(sa.Float)
    conservation_threshold: Mapped[float | None] = mapped_column(
        sa.Float, nullable=True
    )

    host_codon_table: Mapped["CodonTable"] = relationship(
        "CodonTable", foreign_keys=[host_codon_table_id]
    )


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
    identity_percentage: Mapped[float] = mapped_column(sa.Float)
    input_profiles: Mapped[dict] = mapped_column(JSONB, default=lambda: {})
    output_profiles: Mapped[dict] = mapped_column(JSONB, default=lambda: {})
