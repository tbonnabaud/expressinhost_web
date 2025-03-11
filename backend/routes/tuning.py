import asyncio
import re
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
from ..schemas import (
    FineTuningMode,
    ProgressState,
    RunInfoForm,
    RunTuningForm,
    TuningOutput,
)

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


async def stream_sequence_tuning(token: OptionalTokenDependency, form: RunTuningForm):
    # Extract native codon table IDs of the form
    native_codon_table_ids = get_native_codon_table_ids(
        form.nucleotide_file_content, form.sequences_native_codon_tables
    )
    total_sequence_number = len(native_codon_table_ids)

    tuning_state = TuningState().start()
    # Total of steps is the number of sequences
    # plus codon table processing step plus end step
    tuning_state.set_total(total_sequence_number + 2)
    run_start_date = datetime.now(UTC)

    try:
        with context_get_session() as session:
            # Get user if token
            user = get_current_user(session, token) if token else None

            if user is None and isinstance(
                form.five_prime_region_tuning, FineTuningMode
            ):
                form.five_prime_region_tuning = None

            tuning_state.next_step("Process codon tables...")
            yield tuning_state.model_dump_json()

            # Processed codon tables
            native_codon_tables, host_codon_table = get_processed_tables(
                session,
                native_codon_table_ids,
                form.host_codon_table_id,
                form.slow_speed_threshold,
            )

        await asyncio.sleep(0.5)

        sequence_tuner = SequenceTuner(
            form.nucleotide_file_content,
            form.clustal_file_content,
            native_codon_tables,
            host_codon_table,
        )

        tuned_sequences = []
        pipeline = sequence_tuner.tuning_pipeline(
            form.mode, form.conservation_threshold, form.five_prime_region_tuning
        )

        for cpt in range(1, total_sequence_number + 1):
            # Update state before process
            tuning_state.next_step(f"Process sequence {cpt}/{total_sequence_number}...")
            yield tuning_state.model_dump_json()

            try:
                tuned_sequence = next(pipeline)
                tuned_sequences.append(tuned_sequence)
                await asyncio.sleep(0.3)

            except StopIteration:
                break

        tuning_state.next_step("Serialization of the results...")
        yield tuning_state.model_dump_json()

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
                    sequence_number=len(tuned_sequences),
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

        await asyncio.sleep(0.5)
        tuning_state.success()
        tuning_state.result = tuning_output
        yield tuning_state.model_dump_json()

    except (ExpressInHostError, HTTPException) as exc:
        await asyncio.sleep(0.3)
        tuning_state.error(str(exc))
        yield tuning_state.model_dump_json()

    except asyncio.CancelledError:
        await asyncio.sleep(0.3)
        logger.warning("Cancelled.")
        tuning_state.error("Cancelled.")
        yield tuning_state.model_dump_json()

    except Exception as exc:
        await asyncio.sleep(0.3)
        tuning_state.error("Server error.")
        logger.error(exc)
        yield tuning_state.model_dump_json()


@router.post("/run-tuning")
async def run_tuning(
    token: OptionalTokenDependency,
    form: RunTuningForm,
):
    return StreamingResponse(
        stream_sequence_tuning(token, form), media_type="text/event-stream"
    )
