from datetime import datetime
from uuid import UUID

from pydantic import BaseModel, Field


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


class TuningParameters(BaseModel):
    slow_speed_threshold: float
    conservation_threshold: float
    sequence_table_mapping: dict[str, str]


class User(BaseModel):
    id: UUID
    creation_date: datetime
    email: str
    # password: str


class UserForm(BaseModel):
    email: str
    password: str


class CodonTable(BaseModel):
    name: str
    user_id: UUID | None
    creation_date: datetime
    organism: str
    custom: bool


class DefaultCodonTableForm(BaseModel):
    name: str
    user_id: None
    organism: str
    custom: bool


class UserCodonTableForm(BaseModel):
    name: str
    user_id: UUID
    organism: str
    custom: bool


class CodonTranslation(BaseModel):
    codon_table_name: str
    codon: str
    anti_codon: str
    amino_acid: str
    trna_gcn: float
    corresponding_codon: str
    wobble_rate: float


class Result(BaseModel):
    id: UUID
    user_id: UUID | None
    creation_date: datetime
    host_organism: str
    mode: str
    parameters: TuningParameters


class TunedSequence(BaseModel):
    id: UUID
    result_id: UUID
    name: str
    input: str
    output: str
    percentage_identity: float
