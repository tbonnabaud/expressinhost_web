import random

import polars as pl

from .constantes import CONSERVATION_THRESHOLD
from .utils import write_text_to_file

STOP_CODON = "UAA"


def select_conserved_codons(
    potential_codons: list[str],
    potential_ranks: list[float],
    rank: float,
):
    """Find all codons among the potential ones whose rank is at this distance to the native codon."""
    selected_codons = []
    # Find the smallest distance of the potential ranks to the rank of the native codon
    min_delta = min([abs(r - rank) for r in potential_ranks])

    for i in range(len(potential_ranks)):
        # If distance to the Rank is equal the the smallest distance
        if abs(potential_ranks[i] - rank) == min_delta:
            # Subset that contains only the codons that can be corresponding codons (they are "equivalent" codons)
            selected_codons.append(potential_codons[i])

    return selected_codons


def select_codons_for_optimization(
    potential_codons: list[str], potential_ranks: list[float]
) -> list[str]:
    """Find the codons of rank 1 that code for the same amino acid."""
    selected_codons = []

    for i in range(len(potential_ranks)):
        # Subset that contains only the codons that can be corresponding codons. They are "equivalent" codons
        if potential_ranks[i] == 1:
            selected_codons.append(potential_codons[i])

    return selected_codons


def get_random_equivalent_codon(selected_codons: list[str]) -> str:
    """Select randomly the corresponding codon among the equivalent ones."""
    random_index = random.randrange(0, len(selected_codons))

    return selected_codons[random_index]


def direct_mapping(
    cleared_nucleotide_sequences: list[str],
    native_codon_tables: list[pl.DataFrame],
    host_codon_table: pl.DataFrame,
) -> list[str]:
    results = []
    sub_host_table_rows = host_codon_table.select(["AA", "Codon", "Rank"]).rows(
        named=True
    )

    for seq, native_codon_table in zip(
        cleared_nucleotide_sequences, native_codon_tables
    ):
        new_line = ""
        sub_native_rows = native_codon_table.select(["AA", "Codon", "Rank"]).rows(
            named=True
        )

        for t in range(int(len(seq) / 3)):
            input_codon = seq[3 * t : 3 * t + 3]

            if input_codon == "---":
                new_line += input_codon

            else:
                is_codon_in_table = False
                rank = 0.0
                amino_acid = ""

                for native_row in sub_native_rows:
                    # When native codon is found
                    if native_row["Codon"] == input_codon:
                        rank = native_row["Rank"]
                        amino_acid = native_row["AA"]
                        is_codon_in_table = True
                        break

                if not is_codon_in_table:
                    new_line += STOP_CODON

                else:
                    potential_codons = []
                    potential_ranks = []

                    for host_row in sub_host_table_rows:
                        if host_row["AA"] == amino_acid:
                            potential_codons.append(host_row["Codon"])
                            potential_ranks.append(host_row["Rank"])

                    selected_codons = select_conserved_codons(
                        potential_codons, potential_ranks, rank
                    )

                    output_codon = get_random_equivalent_codon(selected_codons)

                    # Add it to the processed sequence
                    new_line += output_codon

        results.append(new_line)

    write_text_to_file("\n".join(results), "tmp/modif_sequences_6.txt")

    return results


def optimisation_and_conservation_1(
    aligned_nucleotide_sequences: list[str],
    symbol_sequence: str,
    native_codon_tables: list[pl.DataFrame],
    host_codon_table: pl.DataFrame,
):
    results = []
    sub_host_table_rows = host_codon_table.select(["AA", "Codon", "Rank"]).rows(
        named=True
    )

    for seq, native_codon_table in zip(
        aligned_nucleotide_sequences, native_codon_tables
    ):
        new_line = ""
        sub_native_rows = native_codon_table.select(["AA", "Codon", "Rank"]).rows(
            named=True
        )

        for t in range(int(len(seq) / 3)):
            input_codon = seq[3 * t : 3 * t + 3]

            if input_codon == "---":
                new_line += input_codon

            else:
                is_codon_in_table = False
                rank = 0.0
                amino_acid = ""

                for native_row in sub_native_rows:
                    # When native codon is found
                    if native_row["Codon"] == input_codon:
                        rank = native_row["Rank"]
                        amino_acid = native_row["AA"]
                        is_codon_in_table = True
                        break

                if not is_codon_in_table:
                    new_line += STOP_CODON

                else:
                    potential_codons = []
                    potential_ranks = []

                    for host_row in sub_host_table_rows:
                        if host_row["AA"] == amino_acid:
                            potential_codons.append(host_row["Codon"])
                            potential_ranks.append(host_row["Rank"])

                    # Check whether the codon has to be optimised or conserved
                    symbol = symbol_sequence[t]

                    # If the symbole for that codon is an asterix we mimic the speed of the native codon
                    if symbol == "*":
                        selected_codons = select_conserved_codons(
                            potential_codons, potential_ranks, rank
                        )

                    # If codon has to be optimized
                    else:
                        selected_codons = select_codons_for_optimization(
                            potential_codons, potential_ranks
                        )

                    output_codon = get_random_equivalent_codon(selected_codons)

                    # Add it to the processed sequence
                    new_line += output_codon

        results.append(new_line)

    write_text_to_file("\n".join(results), "tmp/modif_sequences_6.txt")

    return results


def optimisation_and_conservation_2(
    aligned_nucleotide_sequences: list[str],
    symbol_sequence: str,
    native_codon_tables: list[pl.DataFrame],
    host_codon_table: pl.DataFrame,
):
    cpt_symbols = [0.0 for _ in range(len(aligned_nucleotide_sequences[0]))]

    # Make a loop on all natives to create the corresponding "0" and "S" file ("S" indicate slow codons)
    for seq, native_codon_table in zip(
        aligned_nucleotide_sequences, native_codon_tables
    ):
        sub_native_rows = native_codon_table.select(["Codon", "Symbol_Speed"]).rows(
            named=True
        )

        for t in range(int(len(seq) / 3)):
            input_codon = seq[3 * t : 3 * t + 3]

            if input_codon != "---":
                for native_row in sub_native_rows:
                    # When native codon is found
                    if (
                        native_row["Codon"] == input_codon
                        and native_row["Symbol_Speed"] == "S"
                    ):
                        cpt_symbols[t] += 1.0

    symbol_sequence = ""
    native_table_nb = len(aligned_nucleotide_sequences)

    # Analyse the conservation of slow codons for all native sequences and create a one line sequence
    # of "0" and "S" to indicate where codon speed should be conserved instead of optimised
    for symbol in cpt_symbols:
        if symbol >= CONSERVATION_THRESHOLD * native_table_nb:
            symbol_sequence += "S"
        else:
            symbol_sequence += "0"

    write_text_to_file(symbol_sequence, "tmp/modif_sequences_7.txt")

    # In a similar fashion as in optimisation_and_conservation_1,
    # optimise all sequences but mimic native speed where slow codons are conserved
    results = []
    sub_host_table_rows = host_codon_table.select(["AA", "Codon", "Rank"]).rows(
        named=True
    )

    for seq, native_codon_table in zip(
        aligned_nucleotide_sequences, native_codon_tables
    ):
        new_line = ""
        sub_native_rows = native_codon_table.select(["AA", "Codon", "Rank"]).rows(
            named=True
        )

        for t in range(int(len(seq) / 3)):
            input_codon = seq[3 * t : 3 * t + 3]

            if input_codon == "---":
                new_line += input_codon

            else:
                is_codon_in_table = False
                rank = 0.0
                amino_acid = ""

                for native_row in sub_native_rows:
                    # When native codon is found
                    if native_row["Codon"] == input_codon:
                        rank = native_row["Rank"]
                        amino_acid = native_row["AA"]
                        is_codon_in_table = True
                        break

                if not is_codon_in_table:
                    new_line += STOP_CODON

                else:
                    potential_codons = []
                    potential_ranks = []

                    for host_row in sub_host_table_rows:
                        if host_row["AA"] == amino_acid:
                            potential_codons.append(host_row["Codon"])
                            potential_ranks.append(host_row["Rank"])

                    # Check whether the codon speed has to be optimised or mimicked
                    # If codon has to be optimised
                    if symbol_sequence[t] == "0":
                        selected_codons = select_codons_for_optimization(
                            potential_codons, potential_ranks
                        )

                    else:
                        selected_codons = select_conserved_codons(
                            potential_codons, potential_ranks, rank
                        )

                    output_codon = get_random_equivalent_codon(selected_codons)

                    # Add it to the processed sequence
                    new_line += output_codon

        results.append(new_line)

    write_text_to_file("\n".join(results), "tmp/modif_sequences_6.txt")

    return results
