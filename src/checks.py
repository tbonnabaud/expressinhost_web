from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def check_nucleotides_clustal_identity(
    nucleotide_sequences: list[SeqRecord], clustal_sequences: list[SeqRecord]
) -> bool:
    for nucleotide_record, clustal_record in zip(
        nucleotide_sequences,
        clustal_sequences,
    ):
        try:
            # If cds=True, this checks the sequence starts with a valid alternative start codon
            # (which will be translated as methionine, M), that the sequence length is a multiple
            # of three, and that there is a single in frame stop codon at the end (this will be
            # excluded from the protein sequence, regardless of the to_stop option).
            # If these tests fail, an exception is raised.
            translation = str(nucleotide_record.translate(cds=True).seq)
            clustal_seq = str(clustal_record.seq).replace("-", "")

            if translation != clustal_seq:
                return False

        except TranslationError as exc:
            print(exc)
            return False

    return True


def check_amino_acido_conservation(
    nucleotide_sequences: list[SeqRecord], processed_nucleotide_sequences: list[str]
) -> bool:
    for input_record, processed_seq in zip(
        nucleotide_sequences, processed_nucleotide_sequences
    ):
        input_aa_sequence = input_record.translate(cds=True).seq
        output_aa_sequence = Seq(processed_seq).translate(cds=True)

        if input_aa_sequence != output_aa_sequence:
            print(input_aa_sequence)
            print(output_aa_sequence)
            return False

    return True
