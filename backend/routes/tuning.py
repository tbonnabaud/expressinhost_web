import re
import time
from datetime import UTC, datetime
from uuid import UUID

from fastapi import APIRouter
from fastapi.exceptions import HTTPException
from fastapi.responses import StreamingResponse

from ..authentication import OptionalTokenDependency, get_current_user
from ..core.codon_tables import process_raw_codon_table
from ..core.exceptions import ExpressInHostError
from ..core.sequence_tuning import SequenceTuner
from ..crud.codon_tables import CodonTableRepository
from ..crud.codon_translations import CodonTranslationRepository
from ..crud.results import ResultRepository
from ..crud.run_infos import RunInfoRepository
from ..crud.tuned_sequences import TunedSequenceRepository
from ..database import Session, context_get_session, context_get_session_with_commit
from ..logger import logger
from ..schemas import ProgressState, RunInfoForm, RunTuningForm, TuningOutput

router = APIRouter(tags=["Tuning"])


class TuningState(ProgressState):
    result: TuningOutput | None = None


def get_native_codon_table_ids(
    nucleotide_file_content: str, sequences_native_codon_tables: dict[str, UUID]
) -> list[UUID]:
    """Get the IDs from the mapping of sequence names and ensure the order of FASTA file."""
    sequence_names = re.findall(r"^\>(\S+ ?.*)", nucleotide_file_content, re.MULTILINE)

    return [sequences_native_codon_tables[name] for name in sequence_names]


def process_codon_table_from_db(
    codon_translation_repo: CodonTranslationRepository,
    codon_table_id: UUID,
    slow_speed_threshold: float,
):
    condon_translations = codon_translation_repo.list_from_table(codon_table_id)

    return process_raw_codon_table(
        [tr.__dict__ for tr in condon_translations], slow_speed_threshold
    )


def get_processed_tables(
    session: Session,
    native_codon_table_ids: list[UUID],
    host_codon_table_id: UUID,
    slow_speed_threshold: float,
):
    codon_translation_repo = CodonTranslationRepository(session)
    native_codon_tables = [
        process_codon_table_from_db(
            codon_translation_repo, codon_table_id, slow_speed_threshold
        )
        for codon_table_id in native_codon_table_ids
    ]

    host_codon_table = process_codon_table_from_db(
        codon_translation_repo,
        host_codon_table_id,
        slow_speed_threshold,
    )

    return native_codon_tables, host_codon_table


def stream_sequence_tuning(token: OptionalTokenDependency, form: RunTuningForm):
    tuning_state = TuningState().start()
    tuning_state.set_total(4)
    run_start_date = datetime.now(UTC)

    try:
        with context_get_session() as session:
            # Get user if token
            user = get_current_user(session, token) if token else None
            tuning_state.next_step("Process codon tables...")
            yield tuning_state.model_dump_json()

            # Extract native codon table IDs of the form
            native_codon_table_ids = get_native_codon_table_ids(
                form.nucleotide_file_content, form.sequences_native_codon_tables
            )
            # Processed codon tables
            native_codon_tables, host_codon_table = get_processed_tables(
                session,
                native_codon_table_ids,
                form.host_codon_table_id,
                form.slow_speed_threshold,
            )

        time.sleep(0.5)
        tuning_state.next_step("Preprocess sequences...")
        yield tuning_state.model_dump_json()

        sequence_tuner = SequenceTuner(
            form.nucleotide_file_content,
            form.clustal_file_content,
            native_codon_tables,
            host_codon_table,
        )

        time.sleep(0.5)
        tuning_state.next_step("Process sequences...")
        yield tuning_state.model_dump_json()

        processed_sequences = sequence_tuner.process(
            form.mode, form.conservation_threshold
        )

        time.sleep(0.5)
        tuning_state.next_step("Postprocess sequences...")
        yield tuning_state.model_dump_json()
        tuned_sequences = sequence_tuner.postprocess(processed_sequences)

        # Stringify UUIDs to make the dictionary serializable
        serializable_sequences_native_codon_tables = {
            key: str(value) for key, value in form.sequences_native_codon_tables.items()
        }

        run_end_date = datetime.now(UTC)

        result = {
            "creation_date": run_end_date,
            "name": form.name,
            "host_codon_table_id": form.host_codon_table_id,
            "sequences_native_codon_tables": serializable_sequences_native_codon_tables,
            "mode": form.mode,
            "slow_speed_threshold": form.slow_speed_threshold,
            "conservation_threshold": form.conservation_threshold,
            "five_prime_region_tuning": (
                form.five_prime_region_tuning.model_dump()
                if form.five_prime_region_tuning
                else None
            ),
        }

        with context_get_session_with_commit() as session:
            RunInfoRepository(session).add(
                RunInfoForm(
                    creation_date=run_start_date,
                    duration=(run_end_date - run_start_date),
                    sequence_number=len(processed_sequences),
                    mode=form.mode,
                    slow_speed_threshold=form.slow_speed_threshold,
                    conservation_threshold=form.conservation_threshold,
                    five_prime_region_tuning_mode=(
                        form.five_prime_region_tuning.mode
                        if form.five_prime_region_tuning
                        else None
                    ),
                ).model_dump()
            )

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
            tuning_output = TuningOutput.model_validate(
                {
                    "result": result,
                    "tuned_sequences": tuned_sequences,
                }
            )

        time.sleep(0.5)
        tuning_state.success()
        tuning_state.result = tuning_output
        yield tuning_state.model_dump_json()

    except (ExpressInHostError, HTTPException) as exc:
        time.sleep(0.5)
        tuning_state.error(str(exc))
        yield tuning_state.model_dump_json()

    except Exception as exc:
        time.sleep(0.5)
        tuning_state.error("Server error.")
        logger.error(exc)
        yield tuning_state.model_dump_json()


@router.post("/run-tuning")
def run_tuning(
    token: OptionalTokenDependency,
    form: RunTuningForm,
):
    return StreamingResponse(
        stream_sequence_tuning(token, form), media_type="text/event-stream"
    )
