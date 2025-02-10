from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..logger import logger


def check_nucleotides_clustal_identity(
    nucleotide_record: SeqRecord, clustal_record: SeqRecord
) -> bool:
    """Ensure the translation of nucleotide sequence are matching the amino-acid sequence in the CLUSTAL file."""
    try:
        # Boolean, indicates this is a complete CDS.
        # If cds=True, this checks the sequence starts with a valid alternative start codon
        # (which will be translated as methionine, M), that the sequence length is a multiple
        # of three, and that there is a single in frame stop codon at the end (this will be
        # excluded from the protein sequence, regardless of the to_stop option).
        # If these tests fail, an exception is raised.
        translated_seq = str(nucleotide_record.translate(cds=True).seq)
        clustal_seq = str(clustal_record.seq).replace("-", "")

        return translated_seq == clustal_seq

    except TranslationError as exc:
        logger.error(exc)
        return False


def check_amino_acido_conservation(
    nucleotide_record: SeqRecord,
    processed_nucleotide_sequence: str,
) -> bool:
    input_aa_sequence = nucleotide_record.translate(to_stop=True).seq
    output_aa_sequence = Seq(processed_nucleotide_sequence).translate(to_stop=True)

    return input_aa_sequence == output_aa_sequence
