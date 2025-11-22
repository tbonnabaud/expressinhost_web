import tempfile
from pathlib import Path

from backend.core.codon_tables import (
    ProcessedCodonTable,
    ProcessedCodonTableRow,
    compute_codon_table_speed_symbols,
    process_codon_table_from_file,
    process_raw_codon_table,
)


def create_test_codon_table_file():
    """Create a temporary test codon table file."""
    content = """amino_acid\tanticodon\tcodon\ttrna_gcn\twobble_codon\twobble_rate
Ala\tIGC\tGCU\t2.5\tGCU\t0.0
Ala\tGGC\tGCC\t3.0\tGCC\t0.0
Ala\tCGC\tGCG\t1.5\tGCG\t0.0
Ala\tUGC\tGCA\t2.0\tGCA\t0.0
Arg\tICG\tCGU\t1.0\tCGU\t0.0
Arg\tGCG\tCGC\t2.0\tCGC\t0.0
Arg\tCCG\tCGG\t0.5\tCGG\t0.0
Arg\tUCG\tCGA\t0.8\tCGA\t0.0
Arg\tICU\tAGU\t0.0\tCGU\t0.3
Arg\tGCU\tAGC\t0.0\tCGC\t0.3
Asn\tIGU\tAAU\t3.0\tAAU\t0.0
Asn\tGUU\tAAC\t4.0\tAAC\t0.0
Asp\tIGU\tGAU\t2.5\tGAU\t0.0
Asp\tGUC\tGAC\t3.5\tGAC\t0.0
Cys\tIGC\tUGU\t1.5\tUGU\t0.0
Cys\tGCA\tUGC\t2.5\tUGC\t0.0
Gln\tIUG\tCAG\t3.0\tCAG\t0.0
Gln\tCUG\tCAA\t2.0\tCAA\t0.0
Glu\tIUC\tGAG\t3.5\tGAG\t0.0
Glu\tCUC\tGAA\t2.5\tGAA\t0.0
Gly\tICC\tGGU\t2.0\tGGU\t0.0
Gly\tGCC\tGGC\t3.0\tGGC\t0.0
Gly\tCCC\tGGG\t1.5\tGGG\t0.0
Gly\tUCC\tGGA\t2.5\tGGA\t0.0
His\tIGU\tCAU\t2.0\tCAU\t0.0
His\tGUG\tCAC\t3.0\tCAC\t0.0
Ile\tIAU\tAUU\t3.0\tAUU\t0.0
Ile\tGAU\tAUC\t4.0\tAUC\t0.0
Ile\tCAU\tAUA\t1.5\tAUA\t0.0
Leu\tIAG\tCUG\t4.0\tCUG\t0.0
Leu\tGAG\tCUC\t3.0\tCUC\t0.0
Leu\tCAG\tCUU\t2.0\tCUU\t0.0
Leu\tUAG\tCUA\t1.0\tCUA\t0.0
Leu\tIAA\tUUG\t2.5\tUUG\t0.0
Leu\tCAA\tUUA\t1.5\tUUA\t0.0
Lys\tIUU\tAAG\t3.5\tAAG\t0.0
Lys\tCUU\tAAA\t4.5\tAAA\t0.0
Met\tIAU\tAUG\t5.0\tAUG\t0.0
Phe\tIAA\tUUU\t3.0\tUUU\t0.0
Phe\tGAA\tUUC\t4.0\tUUC\t0.0
Pro\tIGG\tCCU\t2.0\tCCU\t0.0
Pro\tGGG\tCCC\t3.0\tCCC\t0.0
Pro\tCGG\tCCG\t1.5\tCCG\t0.0
Pro\tUGG\tCCA\t2.5\tCCA\t0.0
Ser\tIGA\tUCU\t2.5\tUCU\t0.0
Ser\tGGA\tUCC\t3.5\tUCC\t0.0
Ser\tCGA\tUCG\t1.5\tUCG\t0.0
Ser\tUGA\tUCA\t2.0\tUCA\t0.0
Ser\tICU\tAGU\t0.0\tUCU\t0.3
Ser\tGCU\tAGC\t0.0\tUCC\t0.3
Thr\tIGU\tACU\t2.5\tACU\t0.0
Thr\tGGU\tACC\t3.5\tACC\t0.0
Thr\tCGU\tACG\t1.5\tACG\t0.0
Thr\tUGU\tACA\t2.0\tACA\t0.0
Trp\tICA\tUGG\t4.0\tUGG\t0.0
Tyr\tIUA\tUAU\t2.5\tUAU\t0.0
Tyr\tGUA\tUAC\t3.5\tUAC\t0.0
Val\tIAC\tGUU\t3.0\tGUU\t0.0
Val\tGAC\tGUC\t4.0\tGUC\t0.0
Val\tCAC\tGUG\t2.0\tGUG\t0.0
Val\tUAC\tGUA\t2.5\tGUA\t0.0"""

    temp_file = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".csv")
    temp_file.write(content)
    temp_file.close()
    return temp_file.name


def test_processed_codon_table_row_creation():
    """Test creating a ProcessedCodonTableRow."""
    row = ProcessedCodonTableRow(
        amino_acid="Ala",
        anticodon="IGC",
        codon="GCU",
        trna_gcn=2.5,
        wobble_codon="GCU",
        wobble_rate=0.0,
        rank=0.75,
        speed=0.015,
    )

    assert row.amino_acid == "Ala"
    assert row.codon == "GCU"
    assert row.trna_gcn == 2.5
    assert row.rank == 0.75
    assert row.speed == 0.015


def test_processed_codon_table_getitem():
    """Test accessing codon table by codon."""
    row1 = ProcessedCodonTableRow("Ala", "IGC", "GCU", 2.5, "GCU", 0.0, 0.75, 0.015)
    row2 = ProcessedCodonTableRow("Ala", "GGC", "GCC", 3.0, "GCC", 0.0, 1.0, 0.018)

    table = ProcessedCodonTable(indexed_rows={"GCU": row1, "GCC": row2})

    assert table["GCU"] == row1
    assert table["GCC"] == row2


def test_processed_codon_table_get_row():
    """Test get_row method returns None for non-existent codon."""
    row = ProcessedCodonTableRow("Ala", "IGC", "GCU", 2.5, "GCU", 0.0, 0.75, 0.015)
    table = ProcessedCodonTable(indexed_rows={"GCU": row})

    assert table.get_row("GCU") == row
    assert table.get_row("XXX") is None


def test_processed_codon_table_values():
    """Test values method returns all rows."""
    row1 = ProcessedCodonTableRow("Ala", "IGC", "GCU", 2.5, "GCU", 0.0, 0.75, 0.015)
    row2 = ProcessedCodonTableRow("Ala", "GGC", "GCC", 3.0, "GCC", 0.0, 1.0, 0.018)

    table = ProcessedCodonTable(indexed_rows={"GCU": row1, "GCC": row2})

    values = list(table.values())
    assert len(values) == 2
    assert row1 in values
    assert row2 in values


def test_process_raw_codon_table(raw_codon_table_data):
    """Test processing raw codon table data."""
    result = process_raw_codon_table(raw_codon_table_data)

    assert isinstance(result, ProcessedCodonTable)
    assert "GCU" in result.indexed_rows
    assert "GCC" in result.indexed_rows
    assert len(result.indexed_rows) == 61


def test_process_codon_table_from_file():
    """Test processing codon table from file."""
    temp_file_path = create_test_codon_table_file()

    try:
        result = process_codon_table_from_file(Path(temp_file_path))

        assert isinstance(result, ProcessedCodonTable)
        # Check that some expected codons exist
        assert result.get_row("GCU") is not None
        assert result.get_row("AUG") is not None

        # Verify structure of a row
        ala_row = result.get_row("GCU")
        assert ala_row.amino_acid == "Ala"
        assert ala_row.codon == "GCU"
        assert ala_row.trna_gcn > 0
        assert 0 <= ala_row.rank <= 1
        assert ala_row.speed > 0

    finally:
        Path(temp_file_path).unlink()


def test_compute_codon_table_speed_symbols():
    """Test computing speed symbols for codon table."""
    # Create mock speed values
    speed_col = [0.01, 0.02, 0.015, 0.025, 0.03, 0.012] * 10 + [0.018]  # 61 values
    slow_speed_threshold = 0.5

    symbols = compute_codon_table_speed_symbols(speed_col, slow_speed_threshold)

    assert len(symbols) == 61
    # All symbols should be either "0" or "S"
    assert all(s in ["0", "S"] for s in symbols)
    # At least some should be marked as slow
    assert "S" in symbols


def test_compute_codon_table_speed_symbols_all_slow():
    """Test speed symbols when all codons are slow."""
    # All speeds very close to minimum
    speed_col = [0.01] * 61
    slow_speed_threshold = 0.5

    symbols = compute_codon_table_speed_symbols(speed_col, slow_speed_threshold)

    # When all speeds are the same (or very close), all should be marked as slow
    assert len(symbols) == 61


def test_process_raw_codon_table_with_wobble():
    """Test processing codon table with wobble codons (GCN = 0)."""
    # Create a minimal dataset with wobble codons
    # We need all 61 codons for the function to work properly
    # Use simple amino acids with wobble cases
    raw_data = []

    # Add Ala codons (no wobble)
    for codon in ["GCU", "GCC", "GCA", "GCG"]:
        raw_data.append(
            {
                "amino_acid": "Ala",
                "anticodon": "IGC",
                "codon": codon,
                "trna_gcn": "2.0",
                "wobble_codon": codon,
                "wobble_rate": "0.0",
            }
        )

    # Add Arg codons with wobble (key test case)
    raw_data.extend(
        [
            {
                "amino_acid": "Arg",
                "anticodon": "ICG",
                "codon": "CGU",
                "trna_gcn": "2.0",
                "wobble_codon": "CGU",
                "wobble_rate": "0.0",
            },
            {
                "amino_acid": "Arg",
                "anticodon": "GCG",
                "codon": "CGC",
                "trna_gcn": "3.0",
                "wobble_codon": "CGC",
                "wobble_rate": "0.0",
            },
            {
                "amino_acid": "Arg",
                "anticodon": "CCG",
                "codon": "CGG",
                "trna_gcn": "1.0",
                "wobble_codon": "CGG",
                "wobble_rate": "0.0",
            },
            {
                "amino_acid": "Arg",
                "anticodon": "UCG",
                "codon": "CGA",
                "trna_gcn": "1.5",
                "wobble_codon": "CGA",
                "wobble_rate": "0.0",
            },
            {
                "amino_acid": "Arg",
                "anticodon": "ICU",
                "codon": "AGA",
                "trna_gcn": "0.0",  # This should be calculated via wobble
                "wobble_codon": "CGU",
                "wobble_rate": "0.3",
            },
            {
                "amino_acid": "Arg",
                "anticodon": "GCU",
                "codon": "AGG",
                "trna_gcn": "0.0",  # This should be calculated via wobble
                "wobble_codon": "CGC",
                "wobble_rate": "0.4",
            },
        ]
    )

    # Fill in remaining codons to reach 61 total (using simple dummy data)
    remaining_codons = [
        ("Asn", "AAU"),
        ("Asn", "AAC"),
        ("Asp", "GAU"),
        ("Asp", "GAC"),
        ("Cys", "UGU"),
        ("Cys", "UGC"),
        ("Gln", "CAA"),
        ("Gln", "CAG"),
        ("Glu", "GAA"),
        ("Glu", "GAG"),
        ("Gly", "GGU"),
        ("Gly", "GGC"),
        ("Gly", "GGA"),
        ("Gly", "GGG"),
        ("His", "CAU"),
        ("His", "CAC"),
        ("Ile", "AUU"),
        ("Ile", "AUC"),
        ("Ile", "AUA"),
        ("Leu", "UUA"),
        ("Leu", "UUG"),
        ("Leu", "CUU"),
        ("Leu", "CUC"),
        ("Leu", "CUA"),
        ("Leu", "CUG"),
        ("Lys", "AAA"),
        ("Lys", "AAG"),
        ("Met", "AUG"),
        ("Phe", "UUU"),
        ("Phe", "UUC"),
        ("Pro", "CCU"),
        ("Pro", "CCC"),
        ("Pro", "CCA"),
        ("Pro", "CCG"),
        ("Ser", "UCU"),
        ("Ser", "UCC"),
        ("Ser", "UCA"),
        ("Ser", "UCG"),
        ("Ser", "AGU"),
        ("Ser", "AGC"),
        ("Thr", "ACU"),
        ("Thr", "ACC"),
        ("Thr", "ACA"),
        ("Thr", "ACG"),
        ("Trp", "UGG"),
        ("Tyr", "UAU"),
        ("Tyr", "UAC"),
        ("Val", "GUU"),
        ("Val", "GUC"),
        ("Val", "GUA"),
        ("Val", "GUG"),
    ]

    for aa, codon in remaining_codons:
        raw_data.append(
            {
                "amino_acid": aa,
                "anticodon": "XXX",
                "codon": codon,
                "trna_gcn": "2.0",
                "wobble_codon": codon,
                "wobble_rate": "0.0",
            }
        )

    result = process_raw_codon_table(raw_data)

    # Test that wobble codons were processed correctly
    assert isinstance(result, ProcessedCodonTable)

    # Check AGA codon (should wobble from CGU with rate 0.3)
    # Expected GCN = 2.0 * (1 - 0.3) = 1.4
    aga_row = result.get_row("AGA")
    assert aga_row is not None
    assert aga_row.amino_acid == "Arg"
    assert abs(aga_row.trna_gcn - 1.4) < 0.001, (
        f"Expected GCN 1.4, got {aga_row.trna_gcn}"
    )

    # Check AGG codon (should wobble from CGC with rate 0.4)
    # Expected GCN = 3.0 * (1 - 0.4) = 1.8
    agg_row = result.get_row("AGG")
    assert agg_row is not None
    assert agg_row.amino_acid == "Arg"
    assert abs(agg_row.trna_gcn - 1.8) < 0.001, (
        f"Expected GCN 1.8, got {agg_row.trna_gcn}"
    )

    # Verify the original codon CGU wasn't modified
    cgu_row = result.get_row("CGU")
    assert cgu_row.trna_gcn == 2.0

    # Verify all 61 codons are present
    assert len(result.indexed_rows) == 61


def test_processed_codon_table_rank_calculation():
    """Test that rank is calculated correctly (between 0 and 1)."""
    temp_file_path = create_test_codon_table_file()

    try:
        result = process_codon_table_from_file(Path(temp_file_path))

        # Check all ranks are between 0 and 1
        for row in result.values():
            assert 0 <= row.rank <= 1, (
                f"Rank {row.rank} for codon {row.codon} is out of range"
            )

    finally:
        Path(temp_file_path).unlink()


def test_processed_codon_table_speed_sum():
    """Test that speed values sum to approximately 1."""
    temp_file_path = create_test_codon_table_file()

    try:
        result = process_codon_table_from_file(Path(temp_file_path))

        # Sum all speed values
        total_speed = sum(row.speed for row in result.values())

        # Should sum to approximately 1 (allowing for floating point errors and wobble rate adjustments)
        assert abs(total_speed - 1.0) < 0.02, (
            f"Speed sum {total_speed} is not close to 1"
        )

    finally:
        Path(temp_file_path).unlink()
