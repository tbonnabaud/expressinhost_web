import re
from datetime import UTC, datetime, timedelta
from enum import Enum
from typing import Literal
from uuid import UUID

from pydantic import BaseModel, ConfigDict, Field, computed_field, field_validator


class Status(str, Enum):
    QUEUED = "queued"
    FINISHED = "finished"
    FAILED = "failed"
    STARTED = "started"
    DEFERRED = "deferred"
    SCHEDULED = "scheduled"
    STOPPED = "stopped"
    CANCELED = "canceled"
    NOT_FOUND = "not found"


class ProgressState(BaseModel):
    status: Status
    message: str = ""
    step: int = 0
    total: int = 0
    result: dict | None = None
    exc_info: str | None = None


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


class TuningModeName(str, Enum):
    DIRECT_MAPPING = "direct_mapping"
    OPTIMISATION_AND_CONSERVATION_1 = "optimisation_and_conservation_1"
    OPTIMISATION_AND_CONSERVATION_2 = "optimisation_and_conservation_2"
    PROTEIN_STRUCTURE_ANALYSIS = "protein_structure_analysis"


# Five prime region tuning modes
class PartialUntuningMode(BaseModel):
    mode: Literal["partial_untuning"]
    untuned_codon_number: int


class FineTuningMode(BaseModel):
    """Mode for using OSTIR on five prime region."""

    mode: Literal["fine_tuning"]
    codon_window_size: int
    utr: str

    @field_validator("utr")
    @staticmethod
    def clean(value: str):
        return value.strip()


class SlowedDownMode(BaseModel):
    mode: Literal["slowed_down"]
    slowed_down_codon_number: int


type FivePrimeRegionTuningMode = PartialUntuningMode | FineTuningMode | SlowedDownMode


class RestrictionSite(BaseModel):
    enzyme: str
    sequence: str

    @field_validator("enzyme", "sequence")
    @staticmethod
    def clean(value: str):
        return value.strip()


class RunTuningForm(BaseModel):
    name: str = "Unnamed"
    nucleotide_file_content: str | None
    pdb_file_content: str | None
    clustal_file_content: str | None
    host_codon_table_id: UUID
    sequences_native_codon_tables: dict[str, UUID]
    mode: TuningModeName
    slow_speed_threshold: float | None
    conservation_threshold: float | None
    five_prime_region_tuning: FivePrimeRegionTuningMode | None
    restriction_sites: list[RestrictionSite] | None
    send_email: bool = False

    @field_validator("name")
    @staticmethod
    def clean(value: str):
        return value.strip() if value else "Unnamed"

    @field_validator("nucleotide_file_content", "clustal_file_content")
    @staticmethod
    def remove_carriage_return_characters(value: str | None):
        if value is None:
            return None

        return value.replace("\r", "")


class Result(BaseModel):
    id: UUID | None = None
    user_id: UUID | None = None
    creation_date: datetime = datetime.now(UTC)
    name: str
    host_codon_table_id: UUID
    sequences_native_codon_tables: dict[str, UUID]
    mode: str
    slow_speed_threshold: float | None
    conservation_threshold: float | None
    host_codon_table: CodonTable
    five_prime_region_tuning: FivePrimeRegionTuningMode | None
    # Restriction enzyme recognition sites to avoid
    restriction_sites: list[RestrictionSite] | None


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
    input_profiles: Profiles | None
    output_profiles: Profiles


class TuningOutput(BaseModel):
    result: Result
    tuned_sequences: list[TunedSequence]


class RunInfo(BaseModel):
    id: UUID
    creation_date: datetime
    duration: timedelta
    sequence_number: int
    mode: TuningModeName
    slow_speed_threshold: float | None
    conservation_threshold: float | None
    five_prime_region_tuning_mode: str | None


class RunInfoForm(BaseModel):
    creation_date: datetime
    duration: timedelta
    sequence_number: int
    mode: str
    slow_speed_threshold: float | None
    conservation_threshold: float | None
    five_prime_region_tuning_mode: str | None


class SequenceComparatorForm(BaseModel):
    sequence1: str
    sequence2: str
    host_codon_table_id: UUID
