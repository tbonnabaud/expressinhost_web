import re
from datetime import UTC, datetime
from enum import Enum
from typing import Literal
from uuid import UUID

from pydantic import BaseModel, ConfigDict, Field, computed_field, field_validator


class Status(str, Enum):
    IDLE = "Idle"
    RUNNING = "Running"
    ERROR = "Error"
    SUCCESS = "Success"


class ProgressState(BaseModel):
    status: Status = Status.IDLE
    message: str = ""
    step: int = 0
    total: int = 0

    def start(self, message: str = "Starting..."):
        self.status = Status.RUNNING
        self.step = 0
        self.total = 0
        self.message = message

        return self

    def set_total(self, total: int):
        self.total = total

    def next_step(self, message: str | None = None):
        self.step += 1

        if message is not None:
            self.message = message

    def success(self, message: str = "Success."):
        self.status = Status.SUCCESS
        self.message = message

    def error(self, message: str = "Error."):
        self.status = Status.ERROR
        self.message = message


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
    contact_consent: bool


class UserForm(BaseModel):
    email: str
    password: str
    full_name: str
    contact_consent: bool

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


class UserProfileForm(BaseModel):
    email: str
    full_name: str
    contact_consent: bool

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


class UserPasswordForm(BaseModel):
    reset_token: str
    password: str


class UserRoleForm(BaseModel):
    role: Literal["member", "admin"]


class CodonTranslation(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    # codon_table_id: UUID
    codon: str
    anticodon: str
    amino_acid: str
    trna_gcn: float
    wobble_codon: str
    wobble_rate: float


class CodonTable(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: UUID
    user_id: UUID | None
    creation_date: datetime
    name: str
    organism: str
    source: str | None


class CodonTableWithTranslations(BaseModel):
    id: UUID
    user_id: UUID | None
    creation_date: datetime
    name: str
    organism: str
    source: str | None
    translations: list[CodonTranslation]


class CodonTableFormWithTranslations(BaseModel):
    name: str
    organism: str
    source: Literal["built-in", "user", "Lowe Lab"]
    translations: list[CodonTranslation]

    @field_validator("name", "organism")
    @staticmethod
    def clean(value: str):
        return value.strip()


class UserCodonTableFormWithTranslations(BaseModel):
    name: str
    organism: str
    translations: list[CodonTranslation]

    @field_validator("name", "organism")
    @staticmethod
    def clean(value: str):
        return value.strip()

    @computed_field
    @property
    def source(self) -> str:
        return "user"


class Result(BaseModel):
    id: UUID | None = None
    user_id: UUID | None = None
    creation_date: datetime = datetime.now(UTC)
    name: str
    host_codon_table_id: UUID
    sequences_native_codon_tables: dict[str, UUID]
    mode: str
    slow_speed_threshold: float
    conservation_threshold: float | None
    host_codon_table: CodonTable


class Profiles(BaseModel):
    speed: list[float] | None
    rank: list[float] | None


class TunedSequence(BaseModel):
    id: UUID | None = None
    result_id: UUID | None = None
    name: str
    input: str
    output: str
    identity_percentage: float
    input_profiles: Profiles
    output_profiles: Profiles


class RunTuningForm(BaseModel):
    name: str = "Unnamed"
    nucleotide_file_content: str
    clustal_file_content: str | None
    host_codon_table_id: UUID
    sequences_native_codon_tables: dict[str, UUID]
    mode: Literal[
        "direct_mapping",
        "optimisation_and_conservation_1",
        "optimisation_and_conservation_2",
    ]
    slow_speed_threshold: float
    conservation_threshold: float | None

    @field_validator("name")
    @staticmethod
    def clean(value: str):
        return value.strip() if value else "Unnamed"


class TuningOutput(BaseModel):
    result: Result
    tuned_sequences: list[TunedSequence]
