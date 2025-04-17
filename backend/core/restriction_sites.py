import re

from ..logger import logger
from ..schemas import RestrictionSite
from .codon_tables import ProcessedCodonTable
from .exceptions import ExpressInHostError


def find_recognition_site_positions(
    recognition_site: str, sequence: str
) -> list[tuple[int, int]]:
    """Finds the positions of all occurrences of a recognition site within a given sequence."""
    return [match.span() for match in re.finditer(recognition_site, sequence)]


def replace_codon_by_closest_rank(
    current_codon: str, host_codon_table: ProcessedCodonTable
) -> str:
    """Replace current codon by the one with the closest rank in host codon table."""
    current_codon_row = host_codon_table.get_row(current_codon)
    current_codon_rank = current_codon_row.rank
    # Remove row of current codon and conserve rows with same amino-acids
    filtered_rows = filter(
        lambda x: x.codon != current_codon
        and current_codon_row.amino_acid == x.amino_acid,
        host_codon_table.values(),
    )

    try:
        return min(filtered_rows, key=lambda x: abs(x.rank - current_codon_rank)).codon

    # Avoid the case where only one codon for an amino-acid
    # because it is filtered before and min function requires
    # at least one element
    except ValueError as exc:
        logger.warning(exc)
        return current_codon


def replace_first_codon_within_recognition_site(
    rna_sequence: str,
    recognition_site_positions: list[tuple[int, int]],
    host_codon_table: ProcessedCodonTable,
):
    updated_sequence = []
    treated_site_idx = 0

    for i in range(int(len(rna_sequence) / 3)):
        codon_start = 3 * i
        codon_end = codon_start + 3
        codon = rna_sequence[codon_start:codon_end]

        if treated_site_idx < len(recognition_site_positions):
            site_start, site_end = recognition_site_positions[treated_site_idx]

            if (
                site_start <= codon_start <= site_end
                or site_start + 1
                <= codon_end  # Ensure at least one nucleotide of the codon is inside the site
                <= site_end
            ):
                codon = replace_codon_by_closest_rank(codon, host_codon_table)
                treated_site_idx += 1

        updated_sequence.append(codon)

    return "".join(updated_sequence)


def avoid_restriction_sites_for_sequence(
    rna_sequence: str,
    restriction_sites: list[RestrictionSite],
    host_codon_table: ProcessedCodonTable,
) -> str:
    """Replace codons matching with each enzyme recognition sites."""
    modified_rna_sequence = rna_sequence

    for site in restriction_sites:
        recognition_site_positions = find_recognition_site_positions(
            site.sequence, modified_rna_sequence
        )

        if recognition_site_positions:
            modified_rna_sequence = replace_first_codon_within_recognition_site(
                modified_rna_sequence,
                recognition_site_positions,
                host_codon_table,
            )

    # Check if there are no more enzyme recognition sites
    for site in restriction_sites:
        recognition_site_positions = find_recognition_site_positions(
            site.sequence, modified_rna_sequence
        )

        if recognition_site_positions != []:
            raise ExpressInHostError(
                f"There is still enzyme recognition sites for {site.enzyme} - {site.sequence}."
            )

    return modified_rna_sequence
