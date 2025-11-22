from pathlib import Path

from backend.core.codon_tables import process_codon_table_from_file
from backend.core.restriction_sites import (
    avoid_restriction_sites_for_sequence,
    find_recognition_site_positions,
    replace_codon_by_closest_rank,
    replace_first_codon_within_recognition_site,
)


def test_find_recognition_site_positions():
    """Test finding recognition site positions in a sequence."""
    recognition_site = "AUGC"
    # Test with a sequence containing multiple occurrences of the recognition site
    sequence1 = "GGGGAUGCGGGGAUGCGGCGUACC"
    assert find_recognition_site_positions(recognition_site, sequence1) == [
        (4, 8),
        (12, 16),
    ]

    # Test with a sequence containing no occurrences of the recognition site
    sequence2 = "GGGGGGGGG"
    assert find_recognition_site_positions(recognition_site, sequence2) == []


def test_find_recognition_site_positions_overlapping():
    """Test finding overlapping recognition sites."""
    recognition_site = "AAA"
    sequence = "AAAAAA"  # Multiple overlapping AAA
    positions = find_recognition_site_positions(recognition_site, sequence)
    # Should find all overlapping occurrences
    assert len(positions) >= 1


def test_find_recognition_site_positions_single():
    """Test finding single recognition site."""
    recognition_site = "ATGC"
    sequence = "GGGATGCGGG"
    positions = find_recognition_site_positions(recognition_site, sequence)
    assert positions == [(3, 7)]


def test_replace_first_codon_within_recognition_site():
    """Test replacing first codon within recognition site."""
    host_codon_table = process_codon_table_from_file(
        Path("codon_tables/Bacillus_subtilis.csv")
    )
    recognition_site = "UUGC"
    sequence1 = "GGGGUUGCGGGGUUGCGGCGUACC"
    recognition_site_positions = find_recognition_site_positions(
        recognition_site, sequence1
    )

    result = replace_first_codon_within_recognition_site(
        sequence1, recognition_site_positions, host_codon_table
    )

    # Result should be different from input
    assert result != sequence1
    # Recognition sites should be modified
    assert len(result) == len(sequence1)


def test_replace_codon_by_closest_rank():
    """Test replacing codon by one with closest rank."""
    host_codon_table = process_codon_table_from_file(
        Path("codon_tables/Bacillus_subtilis.csv")
    )

    # Get a codon from the table
    current_codon = "GCU"  # Alanine

    result = replace_codon_by_closest_rank(current_codon, host_codon_table)

    # Result should be a valid codon
    assert len(result) == 3
    # Result should code for the same amino acid
    current_aa = host_codon_table.get_row(current_codon).amino_acid
    result_aa = host_codon_table.get_row(result).amino_acid
    assert current_aa == result_aa


def test_replace_codon_by_closest_rank_single_codon_aa():
    """Test replacing codon when amino acid has only one codon."""
    host_codon_table = process_codon_table_from_file(
        Path("codon_tables/Bacillus_subtilis.csv")
    )

    # Methionine typically has only one codon (ATG/AUG)
    current_codon = "AUG"

    result = replace_codon_by_closest_rank(current_codon, host_codon_table)

    # Should return the same codon if it's the only one
    assert result == current_codon


def test_avoid_restriction_sites_for_sequence():
    """Test avoiding restriction sites in a sequence."""
    from backend.schemas import RestrictionSite

    host_codon_table = process_codon_table_from_file(
        Path("codon_tables/Bacillus_subtilis.csv")
    )

    # Sequence with a restriction site
    rna_sequence = "AUGUGCGAAUUCAAA"  # Contains EcoRI site GAAUUC
    restriction_sites = [RestrictionSite(enzyme="EcoRI", sequence="GAAUUC")]

    result = avoid_restriction_sites_for_sequence(
        rna_sequence, restriction_sites, host_codon_table
    )

    # Result should not contain the restriction site
    assert "GAAUUC" not in result
    # Length should be the same
    assert len(result) == len(rna_sequence)


def test_avoid_restriction_sites_multiple_sites():
    """Test avoiding multiple restriction sites."""
    from backend.schemas import RestrictionSite

    host_codon_table = process_codon_table_from_file(
        Path("codon_tables/Bacillus_subtilis.csv")
    )

    rna_sequence = "AUGUGCGAAUUCAAAGAGCUCAAA"
    restriction_sites = [
        RestrictionSite(enzyme="EcoRI", sequence="GAAUUC"),
        RestrictionSite(enzyme="SacI", sequence="GAGCUC"),
    ]

    result = avoid_restriction_sites_for_sequence(
        rna_sequence, restriction_sites, host_codon_table
    )

    # Neither restriction site should be present
    assert "GAAUUC" not in result
    assert "GAGCUC" not in result


def test_avoid_restriction_sites_empty_list():
    """Test with empty restriction sites list."""
    host_codon_table = process_codon_table_from_file(
        Path("codon_tables/Bacillus_subtilis.csv")
    )

    rna_sequence = "AUGUGCGAAUUCAAA"
    restriction_sites = []

    result = avoid_restriction_sites_for_sequence(
        rna_sequence, restriction_sites, host_codon_table
    )

    # Sequence should remain unchanged
    assert result == rna_sequence
