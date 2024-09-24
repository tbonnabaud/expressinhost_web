import json
import random
from pathlib import Path

import polars as pl

from .checks import check_amino_acido_conservation, check_nucleotides_clustal_identity
from .codon_tables import process_codon_table_from_file
from .exceptions import NoAminoAcidConservation, NoIdenticalSequencesError
from .postprocessing import clear_output_sequences, compare_sequences
from .preprocessing import align_nucleotide_sequences, clear_nucleotide_sequences
from .schemas import TuningParameters
from .utils import (
    find_organism_from_nucleotide_name,
    parse_alignments,
    parse_sequences,
    timeit,
    write_text_to_file,
)

STOP_CODON = "UAA"


def find_amino_acid_and_rank(
    codon: str, table: list[dict]
) -> tuple[str, float] | tuple[None, None]:
    """Return the tuple corresponding to the given codon."""
    for row in table:
        if row["Codon"] == codon:
            return row["AA"], row["Rank"]

    return None, None


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
                amino_acid, rank = find_amino_acid_and_rank(
                    input_codon, sub_native_rows
                )

                if amino_acid is None:
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
) -> list[str]:
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
                amino_acid, rank = find_amino_acid_and_rank(
                    input_codon, sub_native_rows
                )

                if amino_acid is None:
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
    conservation_threshold: float,
) -> list[str]:
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
        if symbol >= conservation_threshold * native_table_nb:
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
                amino_acid, rank = find_amino_acid_and_rank(
                    input_codon, sub_native_rows
                )

                if amino_acid is None:
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


@timeit
def run_tuning(
    nucleotide_file_content: str,
    clustal_file_content: str | None,
    host_organism: str,
    mode: str,
    tuning_parameters: TuningParameters,
) -> tuple[dict[str, str], dict[str, str]]:
    nucleotide_sequences = parse_sequences(nucleotide_file_content, "fasta")
    cleared_nucleotide_sequences = clear_nucleotide_sequences(nucleotide_sequences)

    native_organism_list = [
        find_organism_from_nucleotide_name(record.name)
        for record in nucleotide_sequences
    ]

    native_codon_tables = [
        process_codon_table_from_file(
            name, tuning_parameters.wobble_rate, tuning_parameters.slow_speed_threshold
        )
        for name in native_organism_list
    ]

    host_codon_table = process_codon_table_from_file(
        host_organism,
        tuning_parameters.wobble_rate,
        tuning_parameters.slow_speed_threshold,
    )

    if mode == "direct_mapping":
        output_sequences = direct_mapping(
            cleared_nucleotide_sequences, native_codon_tables, host_codon_table
        )

    else:
        if clustal_file_content is None:
            raise Exception("Clustal file is required.")

        clustal_sequences = parse_sequences(clustal_file_content, "clustal")

        # Ensure sequences are the same in the two files
        match check_nucleotides_clustal_identity(
            nucleotide_sequences, clustal_sequences
        ):
            case (False, errors):
                print(errors)
                raise NoIdenticalSequencesError(
                    "Sequences are not identical. Check their value in the two files and check their order."
                )
            case _:
                pass

        # Only for testing purpose
        write_text_to_file(
            "\n".join([str(r.seq) for r in clustal_sequences]),
            "tmp/modif_sequences_2.txt",
        )

        aligned_nucleotide_sequences = align_nucleotide_sequences(
            clustal_sequences, cleared_nucleotide_sequences
        )

        clustal_alignments = parse_alignments(clustal_file_content, "clustal")[0]
        symbol_sequence = clustal_alignments.column_annotations.get("clustal_consensus")

        if mode == "optimisation_and_conservation_1":
            output_sequences = optimisation_and_conservation_1(
                aligned_nucleotide_sequences,
                symbol_sequence,
                native_codon_tables,
                host_codon_table,
            )

        elif mode == "optimisation_and_conservation_2":
            output_sequences = optimisation_and_conservation_2(
                aligned_nucleotide_sequences,
                symbol_sequence,
                native_codon_tables,
                host_codon_table,
                tuning_parameters.conservation_threshold,
            )

        else:
            raise Exception(
                "Invalid mode. Should be direct_mapping, optimisation_and_conservation_1 or optimisation_and_conservation_2."
            )

    nucleotide_names = [record.name for record in nucleotide_sequences]

    output_path = Path(f"output/{host_organism}")
    output_path.mkdir(parents=True, exist_ok=True)

    cleared_output_sequences = clear_output_sequences(output_sequences)
    write_text_to_file("\n".join(cleared_output_sequences), output_path / f"{mode}.txt")

    identity_percentages = compare_sequences(
        cleared_nucleotide_sequences, cleared_output_sequences
    )

    # Ensure input and ouput nucleotide sequences have the same amino-acid sequences
    match check_amino_acido_conservation(
        nucleotide_sequences, cleared_output_sequences
    ):
        case (False, errors):
            print(errors)
            raise NoAminoAcidConservation(
                "Amino acid sequences from input and output are supposed to be the same."
            )
        case _:
            pass

    output_sequence_mapping = dict(zip(nucleotide_names, cleared_output_sequences))
    identity_percentage_mapping = dict(zip(nucleotide_names, identity_percentages))

    with open(output_path / f"{mode}_identity_percentages.json", "w") as f:
        json.dump(identity_percentage_mapping, f, indent=4)

    return output_sequence_mapping, identity_percentage_mapping
