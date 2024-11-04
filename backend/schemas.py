import re
from datetime import UTC, datetime
from typing import Literal
from uuid import UUID

from pydantic import BaseModel, Field, field_validator


class BadSequence(BaseModel):
    name: str
    msg: str


class MismatchingSequences(BaseModel):
    name: str
    sequence1: str
    sequence2: str


class FilterParams(BaseModel):
    limit: int = Field(100, gt=0, le=100)
    offset: int = Field(0, ge=0)


class Token(BaseModel):
    access_token: str
    token_type: str


class User(BaseModel):
    id: UUID
    creation_date: datetime
    email: str
    role: str
    full_name: str


class UserForm(BaseModel):
    email: str
    password: str
    full_name: str

    @field_validator("email")
    @staticmethod
    def valid_email(value: str):
        value = value.strip().lower()

        if not re.match(r"^[\w\-\.]+@([\w\-]+\.)+\w{2,4}$", value):
            raise ValueError("Invalid e-mail address.")

        return value

    @field_validator("full_name")
    @staticmethod
    def clean(value: str):
        return value.strip()


class CodonTable(BaseModel):
    name: str
    user_id: UUID | None
    creation_date: datetime
    organism: str


class DefaultCodonTableForm(BaseModel):
    name: str
    user_id: None = None
    organism: str

    @field_validator("name", "organism")
    @staticmethod
    def clean(value: str):
        return value.strip()

    @field_validator("organism")
    @staticmethod
    def normalize(value: str):
        # Return value with binomial nomenclature
        return re.sub(r"[\s_\-]+", " ", value).capitalize()


class UserCodonTableForm(BaseModel):
    name: str
    user_id: UUID
    organism: str

    @field_validator("name", "organism")
    @staticmethod
    def clean(value: str):
        return value.strip()

    @field_validator("organism")
    @staticmethod
    def normalize(value: str):
        # Return value with binomial nomenclature
        return re.sub(r"[\s_\-]+", " ", value).capitalize()


class CodonTranslation(BaseModel):
    codon_table_name: str
    codon: str
    anticodon: str
    amino_acid: str
    trna_gcn: float
    corresp_codon: str
    wobble_rate: float


class Result(BaseModel):
    id: UUID | None = None
    user_id: UUID | None = None
    creation_date: datetime = datetime.now(UTC)
    host_codon_table_name: str
    sequences_native_codon_tables: dict[str, str]
    mode: str
    slow_speed_threshold: float
    conservation_threshold: float | None


class TunedSequence(BaseModel):
    id: UUID | None = None
    result_id: UUID | None = None
    name: str
    input: str
    output: str
    identity_percentage: float


class RunTuningForm(BaseModel):
    nucleotide_file_content: str
    clustal_file_content: str | None
    host_codon_table_name: str
    sequences_native_codon_tables: dict[str, str]
    mode: Literal[
        "direct_mapping",
        "optimisation_and_conservation_1",
        "optimisation_and_conservation_2",
    ]
    slow_speed_threshold: float
    conservation_threshold: float | None


class TuningOutput(BaseModel):
    result: Result
    tuned_sequences: list[TunedSequence]
