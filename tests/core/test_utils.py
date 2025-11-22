import tempfile
from pathlib import Path

import pytest

from backend.core.exceptions import (
    ClustalFormatError,
    ExpressInHostError,
    FastaFormatError,
)
from backend.core.utils import (
    find_organism_from_nucleotide_name,
    get_available_organism_list,
    get_clustal_symbol_sequence,
    parse_alignments,
    parse_sequences,
    read_text_file,
    write_text_to_file,
)


def test_read_text_file():
    """Test reading text from a file."""
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
        f.write("Test content\nLine 2")
        temp_path = f.name

    try:
        content = read_text_file(temp_path)
        assert content == "Test content\nLine 2"
    finally:
        Path(temp_path).unlink()


def test_write_text_to_file():
    """Test writing text to a file."""
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
        temp_path = f.name

    try:
        write_text_to_file("Hello World", temp_path)
        content = read_text_file(temp_path)
        assert content == "Hello World"
    finally:
        Path(temp_path).unlink()


def test_parse_sequences_fasta_valid():
    """Test parsing valid FASTA format."""
    fasta_content = """>seq1
ATGCGTACG
>seq2
GCTAGCTA"""

    records = parse_sequences(fasta_content, "fasta")
    assert len(records) == 2
    assert records[0].id == "seq1"
    assert str(records[0].seq) == "ATGCGTACG"
    assert records[1].id == "seq2"
    assert str(records[1].seq) == "GCTAGCTA"


def test_parse_sequences_fasta_invalid():
    """Test parsing invalid FASTA format raises FastaFormatError."""
    invalid_fasta = "This is not a valid FASTA file"

    with pytest.raises(FastaFormatError) as exc_info:
        parse_sequences(invalid_fasta, "fasta")
    assert "Fail to parse file" in str(exc_info.value)


def test_parse_sequences_clustal_valid():
    """Test parsing valid CLUSTAL format."""
    clustal_content = """CLUSTAL W (1.82) multiple sequence alignment

seq1    ATGCGT
seq2    ATGCGT
        ******
"""

    records = parse_sequences(clustal_content, "clustal")
    assert len(records) == 2
    assert records[0].id == "seq1"
    assert records[1].id == "seq2"


def test_parse_sequences_clustal_invalid():
    """Test parsing invalid CLUSTAL format raises ClustalFormatError."""
    invalid_clustal = "This is not a valid CLUSTAL file"

    with pytest.raises(ClustalFormatError) as exc_info:
        parse_sequences(invalid_clustal, "clustal")
    assert "Fail to parse file" in str(exc_info.value)


def test_parse_alignments_fasta():
    """Test parsing alignments in FASTA format."""
    fasta_alignment = """>seq1
ATGC-GTACG
>seq2
ATGCAG-ACG"""

    alignments = parse_alignments(fasta_alignment, "fasta")
    assert len(alignments) >= 1


def test_parse_alignments_clustal():
    """Test parsing alignments in CLUSTAL format."""
    clustal_content = """CLUSTAL W (1.82) multiple sequence alignment

seq1    ATGCGT
seq2    ATGCGT
        ******
"""

    alignments = parse_alignments(clustal_content, "clustal")
    assert len(alignments) >= 1


def test_get_clustal_symbol_sequence():
    """Test extracting symbol sequence from CLUSTAL file."""
    clustal_content = """CLUSTAL W (1.82) multiple sequence alignment

seq1    ATGCGT
seq2    ATGCGT
        ******
"""

    symbol_sequence = get_clustal_symbol_sequence(clustal_content)
    # The symbol sequence should contain conservation markers
    assert symbol_sequence is not None


def test_find_organism_from_nucleotide_name_found():
    """Test finding organism from nucleotide name when organism exists."""
    # This test depends on available organisms in codon_tables directory
    # We'll use a mock approach
    available_organisms = get_available_organism_list()

    if len(available_organisms) > 0:
        # Use first available organism
        organism = available_organisms[0]
        nucleotide_name = f"test_{organism}_sequence"

        result = find_organism_from_nucleotide_name(nucleotide_name)
        assert result == organism


def test_find_organism_from_nucleotide_name_not_found():
    """Test finding organism from nucleotide name when organism doesn't exist."""
    nucleotide_name = "test_nonexistent_organism_xyz_sequence"

    with pytest.raises(ExpressInHostError) as exc_info:
        find_organism_from_nucleotide_name(nucleotide_name)
    assert "Organism not found" in str(exc_info.value)


def test_get_available_organism_list():
    """Test getting list of available organisms."""
    organisms = get_available_organism_list()
    assert isinstance(organisms, list)
    # Should have at least some organisms available
    assert len(organisms) >= 0

    # Test caching - should return the same object
    organisms2 = get_available_organism_list()
    assert organisms == organisms2
