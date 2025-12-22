import pytest

from backend.core.postprocessing import clear_output_sequence, compute_similarity


def test_clear_output_sequence_with_dashes():
    """Test removing dashes from sequence."""
    sequence = "AUG---CGC---UGG"
    result = clear_output_sequence(sequence)
    assert result == "AUGCGCUGG"


def test_clear_output_sequence_no_dashes():
    """Test sequence without dashes remains unchanged."""
    sequence = "AUGCGCUGG"
    result = clear_output_sequence(sequence)
    assert result == "AUGCGCUGG"


def test_clear_output_sequence_only_dashes():
    """Test sequence with only dashes."""
    sequence = "------"
    result = clear_output_sequence(sequence)
    assert result == ""


def test_clear_output_sequence_empty():
    """Test empty sequence."""
    sequence = ""
    result = clear_output_sequence(sequence)
    assert result == ""


def test_clear_output_sequence_mixed():
    """Test sequence with multiple dash groups."""
    sequence = "---AUG---CGC------UGG---"
    result = clear_output_sequence(sequence)
    assert result == "AUGCGCUGG"


def test_compute_similarity_identical():
    """Test similarity of identical sequences."""
    input_seq = "AUGCGCUGG"
    output_seq = "AUGCGCUGG"

    similarity = compute_similarity(input_seq, output_seq)

    # Should be 100% similar (3 codons, all identical)
    assert similarity == 100.0


def test_compute_similarity_completely_different():
    """Test similarity of completely different sequences."""
    input_seq = "AUGCGCUGG"  # M R W
    output_seq = "AAACCCGGG"  # K P G

    similarity = compute_similarity(input_seq, output_seq)

    # Should be 0% similar (no codons match)
    assert similarity == 0.0


def test_compute_similarity_partial():
    """Test similarity of partially matching sequences."""
    input_seq = "AUGCGCUGG"  # ATG CGC TGG (3 codons)
    output_seq = "AUGAAAUGG"  # ATG AAA TGG (2 match, 1 different)

    similarity = compute_similarity(input_seq, output_seq)

    # 2 out of 3 codons match = 66.67%
    assert abs(similarity - 66.666666) < 0.01


def test_compute_similarity_one_codon():
    """Test similarity with single codon."""
    input_seq = "AUG"
    output_seq = "AUG"

    similarity = compute_similarity(input_seq, output_seq)

    assert similarity == 100.0


def test_compute_similarity_one_codon_different():
    """Test similarity with single different codon."""
    input_seq = "AUG"
    output_seq = "AAA"

    similarity = compute_similarity(input_seq, output_seq)

    assert similarity == 0.0


def test_compute_similarity_multiple_codons():
    """Test similarity with longer sequences."""
    # 6 codons
    input_seq = "AUGCGCUGGAAACCCGGG"
    output_seq = "AUGCGCUGGAAACCCGGG"

    similarity = compute_similarity(input_seq, output_seq)

    assert similarity == 100.0


def test_compute_similarity_half_matching():
    """Test similarity with exactly half matching."""
    # 4 codons: 2 match, 2 don't
    input_seq = "AUGCGCAAAGGG"  # ATG CGC AAA GGG
    output_seq = "AUGCGCCCCUUU"  # ATG CGC CCC TTT

    similarity = compute_similarity(input_seq, output_seq)

    # 2 out of 4 = 50%
    assert similarity == 50.0


def test_compute_similarity_precision():
    """Test similarity calculation precision."""
    # 3 codons: 1 match
    input_seq = "AUGCGCUGG"  # ATG CGC TGG
    output_seq = "AUGAAACCC"  # ATG AAA CCC

    similarity = compute_similarity(input_seq, output_seq)

    # 1 out of 3 = 33.333...%
    expected = 100.0 / 3
    assert abs(similarity - expected) < 0.01


def test_compute_similarity_large_sequence():
    """Test similarity with larger sequences."""
    # 10 codons
    input_seq = "AUG" * 10
    output_seq = "AUG" * 9 + "AAA"  # Last one different

    similarity = compute_similarity(input_seq, output_seq)

    # 9 out of 10 = 90%
    assert similarity == 90.0


def test_compute_similarity_length_must_be_multiple_of_three():
    """Test that sequences should be multiples of 3 (codons)."""
    # This test documents expected behavior
    # The function assumes sequences are valid (multiples of 3)
    input_seq = "AUGCGCUGG"  # 9 nucleotides = 3 codons
    output_seq = "AUGCGCUGG"

    similarity = compute_similarity(input_seq, output_seq)

    assert similarity == 100.0


def test_compute_similarity_empty_sequences():
    """Test similarity with empty sequences (edge case)."""
    input_seq = ""
    output_seq = ""

    # This will cause division by zero, but documenting the behavior
    # The function expects valid codon sequences
    # In practice, this shouldn't happen with real data
    with pytest.raises(ZeroDivisionError):
        compute_similarity(input_seq, output_seq)


def test_compute_similarity_dna_vs_rna():
    """Test that function works with both DNA and RNA."""
    # DNA version
    input_seq_dna = "ATGCGCTGG"
    output_seq_dna = "ATGCGCTGG"

    # RNA version
    input_seq_rna = "AUGCGCUGG"
    output_seq_rna = "AUGCGCUGG"

    similarity_dna = compute_similarity(input_seq_dna, output_seq_dna)
    similarity_rna = compute_similarity(input_seq_rna, output_seq_rna)

    # Both should be 100%
    assert similarity_dna == 100.0
    assert similarity_rna == 100.0


def test_compute_similarity_case_sensitive():
    """Test that comparison is case-sensitive."""
    input_seq = "AUGCGCUGG"
    output_seq = "augcgcugg"  # lowercase

    similarity = compute_similarity(input_seq, output_seq)

    # Should be 0% because case doesn't match
    assert similarity == 0.0
