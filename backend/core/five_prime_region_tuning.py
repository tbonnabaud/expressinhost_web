from concurrent.futures import ProcessPoolExecutor
from functools import cache
from itertools import product, repeat

from Bio.Data import CodonTable
from Bio.Seq import Seq
from ostir import run_ostir

from ..logger import logger


@cache
def get_amino_acid_codons_mapping() -> dict[str, list[str]]:
    """Get amino-acid codon mapping from the codon table.

    Returns:
        dict[str, list[str]]: Mapping of amino-acid with its list of codons.
    """
    # importing standard codon table
    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
    # Creating a dictionary mapping amino acids to their codons
    codon_table: dict[str, list[str]] = {}

    for codon, amino_acid in standard_table.forward_table.items():
        if amino_acid not in codon_table:
            codon_table[amino_acid] = []

        codon_table[amino_acid].append(codon)

    return codon_table


def generate_codon_combinations(amino_acid_sequence: str) -> list[str]:
    """Generate all possible codon combinations for a given amino acid sequence."""
    codon_table = get_amino_acid_codons_mapping()

    # Converting the amino acid sequence into a list of codon lists
    codon_lists = [codon_table[amino_acid] for amino_acid in amino_acid_sequence]

    # Generating all possible combinations using itertools.product
    all_codon_combinations = product(*codon_lists)

    # Convert tuples to strings for easier reading
    return ["".join(codons) for codons in all_codon_combinations]


def filter_combinations(
    codon_combinations: list[str], filter_sequences: list[str]
) -> list[str]:
    """Filter out codon combinations that contain any subsequences in the filter list.

    Args:
        codon_combinations (list[str]): List of codon combinations as strings.
        filter_sequences (list[str]): List of subsequences to filter out.

    Returns:
        list[str]: Filtered list of codon combinations.
    """
    return [
        comb
        for comb in codon_combinations
        if not any(seq in comb for seq in filter_sequences)
    ]


def add_UTR_suffix(codon_combinations: list[str], utr: str, suffix: str) -> list[str]:
    """Add UTR and suffix to each codon combination.

    Args:
        codon_combinations (list[str]): List of codon combinations or filtered combinations as strings.
        utr (str): The 5'UTR part of the mRNA entered by user.
        suffix (str): The remaining part of the string (after deduction of codon window).

    Returns:
        list[str]: A new list with the UTR and suffix added to each combination.
    """
    return [f"{utr}{combo}{suffix}" for combo in codon_combinations]


def run_ostir_modified(rna_sequence: str, target_position: int) -> list[dict]:
    """OSTIR identifies all start codon positions and works on all of them.
    Modifying OSTIR output by specifying the start codon at the beginning cleans the output.

    Args:
        rna_sequence (str): mRNA sequence.
        target_position (int): Specifies the start codon position which is after the UTR region.

    Returns:
        list[dict]: List containing dictionaries with 'dG_mRNA'.
    """
    # Simulating the `run_ostir` output
    results = run_ostir(rna_sequence)

    # Return only the first dictionary matching the target position
    for result in results:
        if result["start_position"] == target_position:
            return [result]

    # Return empty list if no match is found
    return []


def find_highest_dG_mRNA(
    results: list[list[dict]], final_combinations: list[str]
) -> str | None:
    """Find the highest 'dG_mRNA' value and its index.
    Use the index to fetch the corresponding sequence from final_combinations.

    Args:
        results (list[list[dict]]): List of lists containing dictionaries with 'dG_mRNA'.
        final_combinations (list[str]): List of sequences.

    Returns:
        (str | None): The sequence from final_combinations corresponding to the highest 'dG_mRNA'.
    """
    max_dG_mRNA = float("-inf")
    max_index = -1

    for i, result in enumerate(results):
        if result:
            dG_mRNA = result[0]["dG_mRNA"]  # Accessing the first dictionary key
            if dG_mRNA > max_dG_mRNA:
                max_dG_mRNA = dG_mRNA
                max_index = i

    if max_index != -1:
        # Get expression if available
        expression_value = results[max_index][0].get("expression", "N/A")
        logger.info(f"Highest dG_mRNA: {max_dG_mRNA}")
        logger.info(f"Expression Value: {expression_value}")
        return final_combinations[max_index]
    else:
        logger.info("No valid dG_mRNA found.")
        return None


def split_open_reading_frame(orf: str, codon_window_size: int) -> tuple[str, str, str]:
    """Splitting the sequence into parts: codon window, suffix, remaining sequence."""
    # Codon window
    codon_window = orf[0 : 3 * codon_window_size]
    # part of the mRNA sequence after the codon window that will be used for OSTIR
    suffix = orf[3 * codon_window_size : 99]
    # remaining part of the mRNA sequence
    remaining_seq = orf[99:]

    return codon_window, suffix, remaining_seq


def optimize_with_ostir(utr: str, orf: str, codon_window_size: int) -> str | None:
    codon_window, suffix, remaining_seq = split_open_reading_frame(
        orf, codon_window_size
    )
    amino_acid_sequence = Seq(codon_window).translate()

    # STEP 1 : Generating Combinations
    codon_combinations = generate_codon_combinations(amino_acid_sequence)

    # STEP 2 : Filter the combinations
    # filter_sequences = ["ACCT", "CCTC", "CTCC", "TCCT"]
    filter_sequences = ["ACCU", "CCUC", "CUCC", "UCCU"]
    filtered_combinations = filter_combinations(codon_combinations, filter_sequences)

    # STEP 3 : Concatenation (Add UTR and suffix to filtered combinations)
    final_combinations = add_UTR_suffix(filtered_combinations, utr, suffix)

    with ProcessPoolExecutor() as executor:
        results = executor.map(
            run_ostir_modified, final_combinations, repeat(len(utr) + 1)
        )
        final_results = [res for res in results if res]

        logger.info(f"The Input sequence: {orf}")
        logger.info(f"Number of combinations: {len(codon_combinations)}")
        logger.info(
            f"Number of combinations after filtering: {len(filtered_combinations)}"
        )

        # Find the sequence with the highest 'dG_mRNA'
        best_sequence = find_highest_dG_mRNA(final_results, final_combinations)
        best_sequence += remaining_seq

        # Remove UTR and return the sequence with the highest 'dG_mRNA'
        return best_sequence[len(utr) :]


# For testing purpose
if __name__ == "__main__":
    ORF = "ATGCGCAGGCTGAAGCGCAAAAGAAGAAAGATGAGGCAGAGGTCCAAGATGGGGAGGGAGAGGAGAGGTGGTAGTGATATGCGCAGGCTGAAGCGCAAAAGAAGAAAGATGAGGCAGAGGTCCAAGATGGGGAGGGAGAGGAGAGGTGGTAGTGATATGCGCAGGCTGAAGCGCAAAAGAAGAAAGATGAGGCAGAGGTCCAAGATGGGGAGGGAGAGGAGAGGTGGTAGTGATATGCGCAGGCTGAAGCGCAAAAGAAGAAAGATGAGGCAGAGGTCCAAGATGGGGAGGGAGAGGAGAGGTGGTAGTGATATGCGCAGGCTGAAGCGCAAAAGAAGAAAGATGAGGCAGAGGTCCAAGATGGGGAGGGAGAGGAGAGGTGGTAGTGATATGCGCAGGCTGAAGCGCAAAAGAAGAAAGATGAGGCAGAGGTCCAAGATGGGGAGGGAGAGGAGAGGTGGTAGTGATATGCGCAGGCTGAAGCGCAAAAGAAGAAAGATGAGGCAGAGGTCCAAGATGGGGAGGGAGAGGAGAGGTGGTAGTGATATGCGCAGGCTGAAGCGCAAAAGAAGAAAGATGAGGCAGAGGTCCAAGATGGGGAGGGAGAGGAGAGGTGGTAGTGATATGCGCAGGCTGAAGCGCAAAAGAAGAAAGATGAGGCAGAGGTCCAAGATGGGGAGGGAGAGGAGAGGTGGTAGTGATATGCGCAGGCTGAAGCGCAAAAGAAGAAAGATGAGGCAGAGGTCCAAGATGGGGAGGGAGAGGAGAGGTGGTAGTGATATGCGCAGGCTGAAGCGCAAAAGAAGAAAGATGAGGCAGAGGTCCAAGATGGGGAGGGAGAGGAGAGGTGGTAGTGATTAA"
    UTR = "ACCCGGCGCTCCATTAAATAGCCGTAGACGGAACTTCGCCTTTCTCTCGGCCTTAGCGCCATTTTTTTGGGTGAGTGTTTTTTGGTTCCTGCGTTGGGATTCCGTGTACAATCCATAGACATCTGACCTCGGCACTTAGCATCATCACAGCAAACTAACTGTAGCCTTTCTCTCTTTCCCTGTAGAAACCTCTGCGCC"

    best_sequence = optimize_with_ostir(UTR, ORF, 5)

    if best_sequence:
        print("Sequence with the highest 'dG_mRNA':", best_sequence)
    else:
        print("No valid sequence found.")
