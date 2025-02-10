import re
from datetime import UTC, datetime
from uuid import UUID

from fastapi import APIRouter

from ..authentication import OptionalTokenDependency, get_current_user
from ..core.codon_tables import process_raw_codon_table
from ..core.sequence_tuning import tune_sequences
from ..crud.codon_tables import CodonTableRepository
from ..crud.codon_translations import CodonTranslationRepository
from ..crud.results import ResultRepository
from ..crud.tuned_sequences import TunedSequenceRepository
from ..database import SessionWithCommitDependency
from ..schemas import RunTuningForm, TuningOutput

router = APIRouter(tags=["Tuning"])


def process_codon_table_from_db(
    codon_translation_repo: CodonTranslationRepository,
    codon_table_id: UUID,
    slow_speed_threshold: float,
):
    condon_translations = codon_translation_repo.list_from_table(codon_table_id)

    return process_raw_codon_table(
        [tr.__dict__ for tr in condon_translations], slow_speed_threshold
    )


@router.post("/run-tuning", response_model=TuningOutput)
def run_tuning(
    session: SessionWithCommitDependency,
    token: OptionalTokenDependency,
    form: RunTuningForm,
):
    sequence_names = re.findall(
        r"^\> *(.*\w)", form.nucleotide_file_content, re.MULTILINE
    )

    codon_translation_repo = CodonTranslationRepository(session)

    native_codon_table_ids = [
        form.sequences_native_codon_tables[name] for name in sequence_names
    ]

    native_codon_tables = [
        process_codon_table_from_db(
            codon_translation_repo, codon_table_id, form.slow_speed_threshold
        )
        for codon_table_id in native_codon_table_ids
    ]

    host_codon_table = process_codon_table_from_db(
        codon_translation_repo,
        form.host_codon_table_id,
        form.slow_speed_threshold,
    )

    tuned_sequences = tune_sequences(
        form.nucleotide_file_content,
        form.clustal_file_content,
        native_codon_tables,
        host_codon_table,
        form.mode,
        form.conservation_threshold,
    )

    # Stringify UUIDs to make the dictionary serializable
    serializable_sequences_native_codon_tables = {
        key: str(value) for key, value in form.sequences_native_codon_tables.items()
    }

    result = {
        "creation_date": datetime.now(UTC),
        "name": form.name,
        "host_codon_table_id": form.host_codon_table_id,
        "sequences_native_codon_tables": serializable_sequences_native_codon_tables,
        "mode": form.mode,
        "slow_speed_threshold": form.slow_speed_threshold,
        "conservation_threshold": form.conservation_threshold,
    }

    user = get_current_user(session, token) if token else None

    if user:
        result["user_id"] = user.id
        result_id = ResultRepository(session).add(result)

        # Attach the tuned sequences to the result
        for seq in tuned_sequences:
            seq["result_id"] = result_id

        TunedSequenceRepository(session).add_batch(tuned_sequences)

    result["host_codon_table"] = CodonTableRepository(session).get(
        user and user.id, form.host_codon_table_id
    )

    return {
        "result": result,
        "tuned_sequences": tuned_sequences,
    }
