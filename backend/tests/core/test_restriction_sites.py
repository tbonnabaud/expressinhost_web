from backend.core.codon_tables import process_codon_table_from_file
from backend.core.restriction_sites import (
    find_recognition_site_positions,
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
    host_codon_table = process_codon_table_from_file(
        "codon_tables/Bacillus_subtilis.csv", 0.5
    )
    recognition_site = "UUGC"
    sequence1 = "GGGGUUGCGGGGUUGCGGCGUACC"
    recognition_site_positions = find_recognition_site_positions(
        recognition_site, sequence1
    )

    assert (
        replace_first_codon_within_recognition_site(
            sequence1, recognition_site_positions, host_codon_table
        )
        == "GGGGUCGCGGGGCUCCGGCGUACC"
    )
