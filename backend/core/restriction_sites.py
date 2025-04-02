import re
from typing import Callable

from .codon_tables import ProcessedCodonTable


def find_recognition_site_positions(
    recognition_site: str, sequence: str
) -> list[tuple[int, int]]:
    """Finds the positions of all occurrences of a recognition site within a given sequence."""
    return [match.span() for match in re.finditer(recognition_site, sequence)]


def replace_codon_by_closest_rank(codon: str, host_codon_table: ProcessedCodonTable):
    "Replace codon by closest rank in host codon table."
    codon_rank = host_codon_table.get_row(codon).rank
    filtered_rows = filter(lambda x: x.codon != codon, host_codon_table.values())

    return min(filtered_rows, key=lambda x: abs(x.rank - codon_rank)).codon


def replace_first_codon_within_recognition_site(
    sequence: str,
    recognition_site_positions: list[tuple[int, int]],
    replacement_strategy: Callable[
        [str], str
    ],  # Function for replacing a codon by another
):
    updated_sequence = []
    treated_site_idx = 0

    for i in range(int(len(sequence) / 3)):
        codon_start = 3 * i
        codon_end = codon_start + 3
        codon = sequence[codon_start:codon_end]

        if treated_site_idx < len(recognition_site_positions):
            site_start, site_end = recognition_site_positions[treated_site_idx]

            if (
                site_start <= codon_start <= site_end
                or site_start + 1 <= codon_end <= site_end
            ):
                codon = replacement_strategy(codon)
                treated_site_idx += 1

        updated_sequence.append(codon)

    return "".join(updated_sequence)
