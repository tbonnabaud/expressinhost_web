from dataclasses import dataclass

from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@dataclass
class BadSequence:
    name: str
    msg: str


@dataclass
class MismatchingSequences:
    name: str
    sequence1: str
    sequence2: str


def check_nucleotides_clustal_identity(
    nucleotide_sequences: list[SeqRecord], clustal_sequences: list[SeqRecord]
) -> tuple[bool, list[MismatchingSequences | BadSequence]]:
    """Ensure the translation of nucleotide sequences are matching the amino-acid sequences in the CLUSTAL file."""
    is_ok = True
    errors = []

    for nucleotide_record, clustal_record in zip(
        nucleotide_sequences,
        clustal_sequences,
    ):
        try:
            # Boolean, indicates this is a complete CDS.
            # If cds=True, this checks the sequence starts with a valid alternative start codon
            # (which will be translated as methionine, M), that the sequence length is a multiple
            # of three, and that there is a single in frame stop codon at the end (this will be
            # excluded from the protein sequence, regardless of the to_stop option).
            # If these tests fail, an exception is raised.
            translated_seq = str(nucleotide_record.translate(cds=True).seq)
            clustal_seq = str(clustal_record.seq).replace("-", "")

            if translated_seq != clustal_seq:
                is_ok = False
                errors.append(
                    MismatchingSequences(
                        nucleotide_record.name, translated_seq, clustal_seq
                    )
                )

        except TranslationError as exc:
            is_ok = False
            errors.append(BadSequence(nucleotide_record.name, str(exc)))

    return is_ok, errors


def check_amino_acido_conservation(
    nucleotide_sequences: list[SeqRecord], processed_nucleotide_sequences: list[str]
) -> tuple[bool, list[MismatchingSequences]]:
    is_ok = True
    errors = []

    for input_record, processed_seq in zip(
        nucleotide_sequences, processed_nucleotide_sequences
    ):
        input_aa_sequence = input_record.translate(to_stop=True).seq
        output_aa_sequence = Seq(processed_seq).translate(to_stop=True)

        if input_aa_sequence != output_aa_sequence:
            is_ok = False
            errors.append(
                MismatchingSequences(
                    input_record.name, str(input_aa_sequence), str(output_aa_sequence)
                )
            )

    return is_ok, errors
