from fastapi import APIRouter

from ..core.postprocessing import compute_similarity
from ..core.sequence_tuning import get_sequence_profiles
from ..core.utils import parse_sequences
from ..crud.codon_tables import CodonTableRepository
from ..database import SessionDependency
from ..schemas import FastaComparatorForm, TunedSequence
from .tuning import process_codon_table_from_db

router = APIRouter(tags=["FASTA comparator"])


@router.get("/compare-fasta", response_model=list[TunedSequence])
def compare_fasta(session: SessionDependency, form: FastaComparatorForm):
    codon_table = process_codon_table_from_db(
        CodonTableRepository(session),
        form.host_codon_table_id,
        form.slow_speed_threshold,
    )

    fasta1_sequences = parse_sequences(form.fasta1, "fasta")
    fasta2_sequences = parse_sequences(form.fasta2, "fasta")

    assert len(fasta1_sequences) == len(fasta2_sequences)

    compared_sequences = []

    for fasta1_seq, fasta2_seq in zip(fasta1_sequences, fasta2_sequences):
        # Stringify sequences
        str_fasta1_seq = str(fasta1_seq.seq)
        str_fasta2_seq = str(fasta2_seq.seq)
        # Computations
        similarity = compute_similarity(str_fasta1_seq, str_fasta2_seq)
        fasta1_seq_profiles = get_sequence_profiles(str_fasta1_seq, codon_table)
        fasta2_seq_profiles = get_sequence_profiles(str_fasta2_seq, codon_table)

        compared_sequences.append(
            TunedSequence(
                name=fasta1_seq.description,
                input=str_fasta1_seq,
                output=str_fasta2_seq,
                similarity=similarity,
                input_profiles=fasta1_seq_profiles,
                output_profiles=fasta2_seq_profiles,
            )
        )

    return compared_sequences
