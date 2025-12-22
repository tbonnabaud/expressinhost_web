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


@router.post("/compare-sequences", response_model=TunedSequence)
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
        CodonTranslationRepository(session), form.host_codon_table_id
    )

    rna_seq1 = form.sequence1.replace("T", "U")
    rna_seq2 = form.sequence2.replace("T", "U")

    assert len(rna_seq1) == len(rna_seq2)

    similarity = compute_similarity(rna_seq1, rna_seq2)
    sequence1_profiles = get_sequence_profiles(rna_seq1, processed_codon_table)
    sequence2_profiles = get_sequence_profiles(rna_seq2, processed_codon_table)

    print(sequence1_profiles)

    return TunedSequence(
        name=f"{codon_table_meta.organism} - {codon_table_meta.name}",
        input=rna_seq1,
        output=rna_seq2,
        identity_percentage=similarity,
        input_profiles=sequence1_profiles,
        output_profiles=sequence2_profiles,
    )
