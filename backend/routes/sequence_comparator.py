from fastapi import APIRouter

from ..authentication import OptionalTokenDependency, get_current_user
from ..core.postprocessing import compute_similarity
from ..core.sequence_tuning import get_sequence_profiles
from ..crud.codon_tables import CodonTableRepository
from ..crud.codon_translations import CodonTranslationRepository
from ..database import SessionDependency
from ..schemas import SequenceComparatorForm, TunedSequence
from .tuning import process_codon_table_from_db

router = APIRouter(tags=["Sequence comparator"])


@router.post("/compare-sequences", response_model=list[TunedSequence])
def compare_sequences(
    session: SessionDependency,
    token: OptionalTokenDependency,
    form: SequenceComparatorForm,
):
    current_user_id = get_current_user(session, token).id if token else None
    codon_table_meta = CodonTableRepository(session).get(
        current_user_id, form.host_codon_table_id
    )
    processed_codon_table = process_codon_table_from_db(
        CodonTranslationRepository(session),
        form.host_codon_table_id,
        form.slow_speed_threshold,
    )

    assert len(form.sequence1) == len(form.sequence2)

    similarity = compute_similarity(form.sequence1, form.sequence2)
    fasta1_seq_profiles = get_sequence_profiles(form.sequence1, processed_codon_table)
    fasta2_seq_profiles = get_sequence_profiles(form.sequence2, processed_codon_table)

    return TunedSequence(
        name=f"{codon_table_meta.organism} - {codon_table_meta.name}",
        input=form.sequence1,
        output=form.sequence2,
        similarity=similarity,
        input_profiles=fasta1_seq_profiles,
        output_profiles=fasta2_seq_profiles,
    )
