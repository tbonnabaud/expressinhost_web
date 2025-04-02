from backend.core.codon_tables import process_codon_table_from_file
from backend.core.restriction_sites import (
    find_recognition_site_positions,
    replace_codon_by_closest_rank,
    replace_first_codon_within_recognition_site,
)


def test_find_recognition_site_positions():
    recognition_site = "AUGC"
    # Uest with a sequence containing multiple occurrences of the recognition site
    sequence1 = "GGGGAUGCGGGGAUGCGGCGUACC"
    assert find_recognition_site_positions(recognition_site, sequence1) == [
        (4, 8),
        (12, 16),
    ]

    # Uest with a sequence containing no occurrences of the recognition site
    sequence2 = "GGGGGGGGG"
    assert find_recognition_site_positions(recognition_site, sequence2) == []


def test_replace_first_codon_within_recognition_site():
    recognition_site = "AUGC"
    sequence1 = "GGGGAUGCGGGGAUGCGGCGUACC"
    recognition_site_positions = find_recognition_site_positions(
        recognition_site, sequence1
    )
    expected_sequence = "GGGXXXGCGGGGXXXCGGCGUACC"

    assert (
        replace_first_codon_within_recognition_site(
            sequence1, recognition_site_positions, lambda _: "XXX"
        )
        == expected_sequence
    )

    host_codon_table = process_codon_table_from_file(
        "codon_tables/Bacillus_subtilis.csv", 0.5
    )

    assert (
        replace_first_codon_within_recognition_site(
            sequence1,
            recognition_site_positions,
            lambda x: replace_codon_by_closest_rank(x, host_codon_table),
        )
        != sequence1
    )
