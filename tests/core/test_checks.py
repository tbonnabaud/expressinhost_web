from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from backend.core.checks import (
    check_amino_acido_conservation,
    check_nucleotides_clustal_identity,
)


def test_check_nucleotides_clustal_identity_matching():
    """Test that matching nucleotide and clustal sequences return True."""
    # Create a nucleotide sequence that translates to "MRW*"
    nucleotide_seq = "ATGCGCTGGTAA"  # ATG=M, CGC=R, TGG=W, TAA=*
    nucleotide_record = SeqRecord(Seq(nucleotide_seq), id="test1")

    # Create matching amino acid sequence
    amino_acid_seq = "MRW"
    clustal_record = SeqRecord(Seq(amino_acid_seq), id="test1")

    result = check_nucleotides_clustal_identity(nucleotide_record, clustal_record)
    assert result is True


def test_check_nucleotides_clustal_identity_with_gaps():
    """Test that sequences with gaps in clustal are handled correctly."""
    # Create a nucleotide sequence
    nucleotide_seq = "ATGCGCTGGTAA"  # ATG=M, CGC=R, TGG=W, TAA=*
    nucleotide_record = SeqRecord(Seq(nucleotide_seq), id="test1")

    # Create amino acid sequence with gaps (should be removed)
    amino_acid_seq = "M-RW"
    clustal_record = SeqRecord(Seq(amino_acid_seq), id="test1")

    result = check_nucleotides_clustal_identity(nucleotide_record, clustal_record)
    assert result is True


def test_check_nucleotides_clustal_identity_not_matching():
    """Test that non-matching sequences return False."""
    # Create a nucleotide sequence
    nucleotide_seq = "ATGCGCTGGTAA"  # ATG=M, CGC=R, TGG=W, TAA=*
    nucleotide_record = SeqRecord(Seq(nucleotide_seq), id="test1")

    # Create different amino acid sequence
    amino_acid_seq = "MKL"
    clustal_record = SeqRecord(Seq(amino_acid_seq), id="test1")

    result = check_nucleotides_clustal_identity(nucleotide_record, clustal_record)
    assert result is False


def test_check_nucleotides_clustal_identity_invalid_cds():
    """Test that invalid CDS (not multiple of 3, no start codon, etc.) returns False."""
    # Create invalid nucleotide sequence (not multiple of 3)
    nucleotide_seq = "ATGCG"
    nucleotide_record = SeqRecord(Seq(nucleotide_seq), id="test1")

    amino_acid_seq = "M"
    clustal_record = SeqRecord(Seq(amino_acid_seq), id="test1")

    result = check_nucleotides_clustal_identity(nucleotide_record, clustal_record)
    assert result is False


def test_check_nucleotides_clustal_identity_no_stop_codon():
    """Test that sequence without stop codon returns False."""
    # Create a nucleotide sequence without stop codon
    nucleotide_seq = "ATGCGCTGG"  # ATG=M, CGC=R, TGG=W (no stop)
    nucleotide_record = SeqRecord(Seq(nucleotide_seq), id="test1")

    amino_acid_seq = "MRW"
    clustal_record = SeqRecord(Seq(amino_acid_seq), id="test1")

    result = check_nucleotides_clustal_identity(nucleotide_record, clustal_record)
    assert result is False


def test_check_amino_acido_conservation_matching():
    """Test that sequences with conserved amino acids return True."""
    # Create a nucleotide sequence
    nucleotide_seq = "ATGCGCTGGTAA"  # ATG=M, CGC=R, TGG=W, TAA=*
    nucleotide_record = SeqRecord(Seq(nucleotide_seq), id="test1")

    # Create processed nucleotide sequence with different codons but same amino acids
    # ATG=M, AGA=R, TGG=W, TAA=*
    processed_nucleotide = "ATGAGATGGTAA"

    result = check_amino_acido_conservation(nucleotide_record, processed_nucleotide)
    assert result is True


def test_check_amino_acido_conservation_not_matching():
    """Test that sequences with different amino acids return False."""
    # Create a nucleotide sequence
    nucleotide_seq = "ATGCGCTGGTAA"  # ATG=M, CGC=R, TGG=W, TAA=*
    nucleotide_record = SeqRecord(Seq(nucleotide_seq), id="test1")

    # Create processed nucleotide sequence with different amino acids
    # ATG=M, AAA=K, CTG=L, TAA=*
    processed_nucleotide = "ATGAAACTGTAA"

    result = check_amino_acido_conservation(nucleotide_record, processed_nucleotide)
    assert result is False


def test_check_amino_acido_conservation_with_rna():
    """Test conservation check with RNA sequence (U instead of T)."""
    # Create a DNA nucleotide sequence
    nucleotide_seq = "ATGCGCTGGTAA"
    nucleotide_record = SeqRecord(Seq(nucleotide_seq), id="test1")

    # Create processed RNA sequence with same amino acids
    processed_nucleotide = "AUGAGAUGGUAA"

    result = check_amino_acido_conservation(nucleotide_record, processed_nucleotide)
    assert result is True


def test_check_amino_acido_conservation_partial_sequence():
    """Test conservation check stops at first stop codon."""
    # Create a nucleotide sequence with stop codon
    nucleotide_seq = "ATGCGCTGGTAA"  # M R W *
    nucleotide_record = SeqRecord(Seq(nucleotide_seq), id="test1")

    # Processed sequence with additional codons after stop (should be ignored)
    processed_nucleotide = "AUGAGAUGGUAAAAAAAG"

    result = check_amino_acido_conservation(nucleotide_record, processed_nucleotide)
    assert result is True
