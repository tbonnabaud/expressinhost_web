import re
import time
from datetime import UTC, datetime
from uuid import UUID

from fastapi import APIRouter, Depends, Request
from fastapi.exceptions import HTTPException
from fastapi.responses import StreamingResponse
from rq import get_current_job

from ..authentication import OptionalTokenDependency, check_is_member, get_current_user
from ..core.codon_tables import ProcessedCodonTable, process_raw_codon_table
from ..core.exceptions import ExpressInHostError
from ..core.sequence_tuning import SequenceTuner, StructureTuner
from ..crud.codon_tables import CodonTableRepository
from ..crud.codon_translations import CodonTranslationRepository
from ..crud.results import ResultRepository
from ..crud.run_infos import RunInfoRepository
from ..crud.tuned_sequences import TunedSequenceRepository
from ..database import Session, context_get_session, context_get_session_with_commit
from ..email_service import send_email
from ..job_manager import (
    cancel_job,
    heavy_queue,
    light_queue,
    stream_job_state,
    update_job_meta,
)
from ..logger import logger
from ..schemas import (
    CodonTable,
    FineTuningMode,
    ProgressState,
    RunInfoForm,
    RunTuningForm,
    TuningModeName,
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
    codon_translation_repo: CodonTranslationRepository, codon_table_id: UUID
) -> ProcessedCodonTable:
    condon_translations = codon_translation_repo.list_from_table(codon_table_id)

    return process_raw_codon_table([tr.__dict__ for tr in condon_translations])


def get_processed_tables(
    session: Session,
    native_codon_table_ids: list[UUID],
    host_codon_table_id: UUID,
):
    codon_translation_repo = CodonTranslationRepository(session)
    native_codon_tables = [
        process_codon_table_from_db(codon_translation_repo, codon_table_id)
        for codon_table_id in native_codon_table_ids
    ]

    host_codon_table = process_codon_table_from_db(
        codon_translation_repo, host_codon_table_id
    )

    return native_codon_tables, host_codon_table


def tune_sequences(token: OptionalTokenDependency, base_url: str, form: RunTuningForm):
    job = get_current_job()

    # Extract native codon table IDs of the form
    native_codon_table_ids = get_native_codon_table_ids(
        form.nucleotide_file_content, form.sequences_native_codon_tables
    )

    if form.mode == TuningModeName.PROTEIN_STRUCTURE_ANALYSIS:
        # Only one sequence of amino-acids, provided by the PDB file
        total_sequence_number = 1
    else:
        total_sequence_number = len(native_codon_table_ids)

    step = 0
    # Total of steps is the number of sequences
    # plus codon table processing step plus end step
    update_job_meta(job, "Init...", step, total_sequence_number + 2)
    run_start_date = datetime.now(UTC)

    try:
        with context_get_session() as session:
            # Get user if token
            user = get_current_user(session, token) if token else None

            if user is None and isinstance(
                form.five_prime_region_tuning, FineTuningMode
            ):
                form.five_prime_region_tuning = None

            step += 1
            update_job_meta(job, "Process codon tables...", step)

            # Processed codon tables
            native_codon_tables, host_codon_table = get_processed_tables(
                session, native_codon_table_ids, form.host_codon_table_id
            )

        time.sleep(0.5)

        tuned_sequences = []

        if form.mode == TuningModeName.PROTEIN_STRUCTURE_ANALYSIS:
            step += 1
            update_job_meta(job, "Process protein structure...", step)
            structure_tuner = StructureTuner(form.pdb_file_content, host_codon_table)
            tuned_sequence = structure_tuner.tuning(
                form.five_prime_region_tuning,
                form.restriction_sites,
            )
            tuned_sequences.append(tuned_sequence)

        else:
            sequence_tuner = SequenceTuner(
                form.nucleotide_file_content,
                form.clustal_file_content,
                native_codon_tables,
                host_codon_table,
            )

            pipeline = sequence_tuner.tuning_pipeline(
                form.mode,
                form.slow_speed_threshold,
                form.conservation_threshold,
                form.five_prime_region_tuning,
                form.restriction_sites,
            )

            for seq_no in range(1, total_sequence_number + 1):
                # Update state before process
                step += 1
                update_job_meta(
                    job, f"Process sequence {seq_no}/{total_sequence_number}...", step
                )

                try:
                    tuned_sequence = next(pipeline)
                    tuned_sequences.append(tuned_sequence)
                    time.sleep(0.3)

                except StopIteration:
                    break

        step += 1
        update_job_meta(job, "Serialization of the results...", step)

        # Stringify UUIDs to make the dictionary serializable
        serializable_sequences_native_codon_tables = {
            key: str(value) for key, value in form.sequences_native_codon_tables.items()
        }

        run_end_date = datetime.now(UTC)

        logger.info(form.restriction_sites)

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
            "restriction_sites": (
                [site.model_dump() for site in form.restriction_sites]
                if form.restriction_sites
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

            host_codon_table = CodonTableRepository(session).get(
                user and user.id, form.host_codon_table_id
            )
            result["host_codon_table"] = CodonTable.model_validate(
                host_codon_table
            ).model_dump()

        time.sleep(0.5)

        if user and form.send_email and result_id:
            result_url = f"{base_url}results/{result_id}"
            send_email(
                user.email,
                f'Job "{form.name}" completed',
                f'Job "{form.name}" completed. Here is the link: {result_url}',
            )

        return {
            "result": result,
            "tuned_sequences": tuned_sequences,
        }

    except (ExpressInHostError, HTTPException) as exc:
        update_job_meta(job, str(exc), step)

        if user and form.send_email:
            send_email(
                user.email,
                f'Job "{form.name}" failed',
                f'Job "{form.name}" failed with the following error message:\n"{exc}"',
            )

        raise exc

    except Exception as exc:
        logger.error(exc)
        update_job_meta(job, "Server error.", step)

        if user and form.send_email:
            send_email(
                user.email,
                f'Job "{form.name}" failed',
                f'Job "{form.name}" failed with the following error message:\n"Server error"',
            )

        raise Exception("Server error.")


@router.post("/run-tuning")
async def run_tuning(
    request: Request,
    token: OptionalTokenDependency,
    form: RunTuningForm,
):
    selected_queue = (
        heavy_queue
        if isinstance(form.five_prime_region_tuning, FineTuningMode)
        else light_queue
    )
    job = selected_queue.enqueue(
        tune_sequences,
        token,
        str(request.base_url),
        form,
        job_timeout=-1,  # -1 for infinite timeout, if None the default value 180 is used
        failure_ttl=500,
    )

    return job.id


@router.get("/tuning/state/{job_id}")
def stream_tuning_state(job_id: str):
    return StreamingResponse(stream_job_state(job_id), media_type="text/event-stream")


@router.delete("/tuning/{job_id}", dependencies=[Depends(check_is_member)])
def cancel_tuning(job_id: str):
    cancel_job(job_id)

    return f"Job {job_id} cancelled."
