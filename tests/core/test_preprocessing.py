from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from backend.core.preprocessing import align_nucleotide_sequences, dna_to_rna_sequences


def test_dna_to_rna_sequences():
    """Test conversion of DNA sequences to RNA."""
    # Create DNA sequence records
    dna_seq1 = SeqRecord(Seq("ATGCGT"), id="seq1")
    dna_seq2 = SeqRecord(Seq("GCTAGC"), id="seq2")

    nucleotide_records = [dna_seq1, dna_seq2]

    result = dna_to_rna_sequences(nucleotide_records)

    assert len(result) == 2
    assert result[0] == "AUGCGU"  # T -> U
    assert result[1] == "GCUAGC"  # T -> U


def test_dna_to_rna_sequences_single():
    """Test conversion of single DNA sequence to RNA."""
    dna_seq = SeqRecord(Seq("TTTTAAAA"), id="seq1")

    result = dna_to_rna_sequences([dna_seq])

    assert len(result) == 1
    assert result[0] == "UUUUAAAA"


def test_dna_to_rna_sequences_empty():
    """Test conversion of empty list."""
    result = dna_to_rna_sequences([])

    assert len(result) == 0
    assert result == []


def test_align_nucleotide_sequences_no_gaps():
    """Test aligning nucleotide sequences when clustal has no gaps."""
    # Create clustal record without gaps
    clustal_record = SeqRecord(Seq("MRW"), id="seq1")  # 3 amino acids
    # Corresponding nucleotide sequence (9 nucleotides + stop codon)
    nucleo_seq = "AUGCGCUGGUAA"

    clustal_records = [clustal_record]
    cleared_nucleotide_sequences = [nucleo_seq]

    result = align_nucleotide_sequences(clustal_records, cleared_nucleotide_sequences)

    assert len(result) == 1
    # Should have the original sequence plus stop codon preserved
    assert result[0] == "AUGCGCUGGUAA"


def test_align_nucleotide_sequences_with_gaps():
    """Test aligning nucleotide sequences when clustal has gaps."""
    # Create clustal record with gaps (- represents gap)
    clustal_record = SeqRecord(Seq("M-RW"), id="seq1")  # Gap after first amino acid
    # Corresponding nucleotide sequence (without gap representation)
    nucleo_seq = "AUGCGCUGGUAA"  # MRW* (9 nt + 3 stop)

    clustal_records = [clustal_record]
    cleared_nucleotide_sequences = [nucleo_seq]

    result = align_nucleotide_sequences(clustal_records, cleared_nucleotide_sequences)

    assert len(result) == 1
    # Gap in amino acid should translate to "---" (3 nucleotides)
    assert "---" in result[0]
    # Should be: ATG(M) + ---(gap) + CGC(R) + TGG(W) + TAA(stop)
    assert result[0] == "AUG---CGCUGGUAA"


def test_align_nucleotide_sequences_multiple_gaps():
    """Test aligning with multiple gaps."""
    # Clustal with multiple gaps
    clustal_record = SeqRecord(Seq("M--R"), id="seq1")
    # Original sequence without gaps: MR* = AUG CGC UAA
    nucleo_seq = "AUGCGCUAA"

    clustal_records = [clustal_record]
    cleared_nucleotide_sequences = [nucleo_seq]

    result = align_nucleotide_sequences(clustal_records, cleared_nucleotide_sequences)

    assert len(result) == 1
    # Should be: ATG(M) + ------(2 gaps) + CGC(R) + TAA(stop)
    assert result[0] == "AUG------CGCUAA"


def test_align_nucleotide_sequences_multiple_sequences():
    """Test aligning multiple sequences."""
    # Two clustal records
    clustal_record1 = SeqRecord(Seq("MR"), id="seq1")
    clustal_record2 = SeqRecord(Seq("M-R"), id="seq2")  # With gap

    # Corresponding nucleotide sequences
    nucleo_seq1 = "AUGCGCUAA"  # MR*
    nucleo_seq2 = "AUGCGCUAA"  # MR* (gap added during alignment)

    clustal_records = [clustal_record1, clustal_record2]
    cleared_nucleotide_sequences = [nucleo_seq1, nucleo_seq2]

    result = align_nucleotide_sequences(clustal_records, cleared_nucleotide_sequences)

    assert len(result) == 2
    assert result[0] == "AUGCGCUAA"  # No gaps
    assert result[1] == "AUG---CGCUAA"  # One gap


def test_align_nucleotide_sequences_gap_at_start():
    """Test alignment with gap at the beginning."""
    clustal_record = SeqRecord(Seq("-MR"), id="seq1")
    nucleo_seq = "AUGCGCUAA"

    clustal_records = [clustal_record]
    cleared_nucleotide_sequences = [nucleo_seq]

    result = align_nucleotide_sequences(clustal_records, cleared_nucleotide_sequences)

    assert len(result) == 1
    assert result[0] == "---AUGCGCUAA"


def test_align_nucleotide_sequences_gap_at_end():
    """Test alignment with gap at the end."""
    clustal_record = SeqRecord(Seq("MR-"), id="seq1")
    nucleo_seq = "AUGCGCUAA"

    clustal_records = [clustal_record]
    cleared_nucleotide_sequences = [nucleo_seq]

    result = align_nucleotide_sequences(clustal_records, cleared_nucleotide_sequences)

    assert len(result) == 1
    # Gap at end, plus original stop codon
    assert result[0] == "AUGCGC---UAA"


def test_align_nucleotide_sequences_preserves_stop_codon():
    """Test that stop codon is preserved at the end."""
    clustal_record = SeqRecord(Seq("MRW"), id="seq1")
    nucleo_seq = "AUGCGCUGGUAG"  # Stop codon UAG

    clustal_records = [clustal_record]
    cleared_nucleotide_sequences = [nucleo_seq]

    result = align_nucleotide_sequences(clustal_records, cleared_nucleotide_sequences)

    assert len(result) == 1
    # Last 3 characters should be the stop codon
    assert result[0].endswith("UAG")
    assert result[0] == "AUGCGCUGGUAG"


def test_align_nucleotide_sequences_empty():
    """Test aligning empty lists."""
    result = align_nucleotide_sequences([], [])

    assert len(result) == 0
    assert result == []
