import re
from typing import Annotated

from fastapi import APIRouter, Depends

from ..authentication import get_current_user, optional_oauth2_scheme
from ..core.codon_tables import process_codon_table_from_file
from ..core.sequence_tuning import run_tuning
from ..core.utils import find_organism_from_nucleotide_name
from ..crud.results import ResultRepository
from ..crud.tuned_sequences import TunedSequenceRepository
from ..database import SessionDependency

# from ..models import CodonTranslation, Result
from ..schemas import RunTuningForm, TuningOutput

router = APIRouter(tags=["Tuning"])


@router.post("/run-tuning", response_model=TuningOutput)
def launch_tuning(
    session: SessionDependency,
    token: Annotated[str | None, Depends(optional_oauth2_scheme)],
    form: RunTuningForm,
):
    sequence_names = re.findall(
        r"^\> *(.*\w)", form.nucleotide_file_content, re.MULTILINE
    )

    native_codon_table_names = [
        find_organism_from_nucleotide_name(name) for name in sequence_names
    ]

    native_codon_tables = [
        process_codon_table_from_file(name, form.slow_speed_threshold)
        for name in native_codon_table_names
    ]

    host_codon_table = process_codon_table_from_file(
        form.host_codon_table_name,
        form.slow_speed_threshold,
    )

    tuned_sequences = run_tuning(
        form.nucleotide_file_content,
        form.clustal_file_content,
        native_codon_tables,
        host_codon_table,
        form.mode,
        form.conservation_threshold,
    )

    result = {
        "host_codon_table_name": form.host_codon_table_name,
        "sequences_native_codon_tables": form.sequences_native_codon_tables,
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

    return {
        "result": result,
        "tuned_sequences": tuned_sequences,
    }
