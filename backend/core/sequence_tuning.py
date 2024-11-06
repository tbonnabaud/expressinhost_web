import random

from .checks import check_amino_acido_conservation, check_nucleotides_clustal_identity
from .codon_tables import ProcessedCodonTable
from .exceptions import NoAminoAcidConservation, NoIdenticalSequencesError
from .postprocessing import clear_output_sequences, compare_sequences
from .preprocessing import align_nucleotide_sequences, clear_nucleotide_sequences
from .utils import (  # write_text_to_file,
    get_clustal_symbol_sequence,
    parse_sequences,
    timeit,
)

STOP_CODON = "UAA"


def find_amino_acid_and_rank(
    codon: str, table: ProcessedCodonTable
) -> tuple[str, float] | tuple[None, None]:
    """Return the tuple corresponding to the given codon."""
    row = table.get(codon)

    if row:
        return row.amino_acid, row.rank

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
    native_codon_tables: list[ProcessedCodonTable],
    host_codon_table: ProcessedCodonTable,
) -> list[str]:
    results = []

    for seq, native_codon_table in zip(
        cleared_nucleotide_sequences, native_codon_tables
    ):
        new_line = ""

        for t in range(int(len(seq) / 3)):
            input_codon = seq[3 * t : 3 * t + 3]

            if input_codon == "---":
                new_line += input_codon

            else:
                # Find amino-acid corresponding to native codon
                amino_acid, rank = find_amino_acid_and_rank(
                    input_codon, native_codon_table
                )

                if amino_acid is None:
                    new_line += STOP_CODON

                else:
                    potential_codons = []
                    potential_ranks = []

                    for host_row in host_codon_table.values():
                        if host_row.amino_acid == amino_acid:
                            potential_codons.append(host_row.codon)
                            potential_ranks.append(host_row.rank)

                    selected_codons = select_conserved_codons(
                        potential_codons, potential_ranks, rank
                    )

                    output_codon = get_random_equivalent_codon(selected_codons)

                    # Add it to the processed sequence
                    new_line += output_codon

        results.append(new_line)

    # write_text_to_file("\n".join(results), "tmp/modif_sequences_6.txt")

    return results


def optimisation_and_conservation_1(
    aligned_nucleotide_sequences: list[str],
    symbol_sequence: str,
    native_codon_tables: list[ProcessedCodonTable],
    host_codon_table: ProcessedCodonTable,
) -> list[str]:
    results = []

    for seq, native_codon_table in zip(
        aligned_nucleotide_sequences, native_codon_tables
    ):
        new_line = ""

        for t in range(int(len(seq) / 3)):
            input_codon = seq[3 * t : 3 * t + 3]

            if input_codon == "---":
                new_line += input_codon

            else:
                # Find amino-acid corresponding to native codon
                amino_acid, rank = find_amino_acid_and_rank(
                    input_codon, native_codon_table
                )

                if amino_acid is None:
                    new_line += STOP_CODON

                else:
                    potential_codons = []
                    potential_ranks = []

                    for host_row in host_codon_table.values():
                        if host_row.amino_acid == amino_acid:
                            potential_codons.append(host_row.codon)
                            potential_ranks.append(host_row.rank)

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

    # write_text_to_file("\n".join(results), "tmp/modif_sequences_6.txt")

    return results


def optimisation_and_conservation_2(
    aligned_nucleotide_sequences: list[str],
    symbol_sequence: str,
    native_codon_tables: list[ProcessedCodonTable],
    host_codon_table: ProcessedCodonTable,
    conservation_threshold: float,
) -> list[str]:
    cpt_symbols = [0.0 for _ in range(len(aligned_nucleotide_sequences[0]))]

    # Make a loop on all natives to create the corresponding "0" and "S" file ("S" indicate slow codons)
    for seq, native_codon_table in zip(
        aligned_nucleotide_sequences, native_codon_tables
    ):
        for t in range(int(len(seq) / 3)):
            input_codon = seq[3 * t : 3 * t + 3]

            if input_codon != "---":
                for native_row in native_codon_table.values():
                    # When native codon is found
                    if (
                        native_row.codon == input_codon
                        and native_row.symbol_speed == "S"
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

    # write_text_to_file(symbol_sequence, "tmp/modif_sequences_7.txt")

    # In a similar fashion as in optimisation_and_conservation_1,
    # optimise all sequences but mimic native speed where slow codons are conserved
    results = []

    for seq, native_codon_table in zip(
        aligned_nucleotide_sequences, native_codon_tables
    ):
        new_line = ""

        for t in range(int(len(seq) / 3)):
            input_codon = seq[3 * t : 3 * t + 3]

            if input_codon == "---":
                new_line += input_codon

            else:
                # Find amino-acid corresponding to native codon
                amino_acid, rank = find_amino_acid_and_rank(
                    input_codon, native_codon_table
                )

                if amino_acid is None:
                    new_line += STOP_CODON

                else:
                    potential_codons = []
                    potential_ranks = []

                    for host_row in host_codon_table.values():
                        if host_row.amino_acid == amino_acid:
                            potential_codons.append(host_row.codon)
                            potential_ranks.append(host_row.rank)

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

    # write_text_to_file("\n".join(results), "tmp/modif_sequences_6.txt")

    return results


@timeit
def run_tuning(
    nucleotide_file_content: str,
    clustal_file_content: str | None,
    native_codon_tables: list[ProcessedCodonTable],
    host_codon_table: ProcessedCodonTable,
    mode: str,
    conservation_threshold: float | None,
) -> list[dict]:
    nucleotide_sequences = parse_sequences(nucleotide_file_content, "fasta")
    cleared_nucleotide_sequences = clear_nucleotide_sequences(nucleotide_sequences)

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

        # write_text_to_file(
        #     "\n".join([str(r.seq) for r in clustal_sequences]),
        #     "tmp/modif_sequences_2.txt",
        # )

        aligned_nucleotide_sequences = align_nucleotide_sequences(
            clustal_sequences, cleared_nucleotide_sequences
        )

        symbol_sequence = get_clustal_symbol_sequence(clustal_file_content)

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
                conservation_threshold,
            )

        else:
            raise Exception(
                "Invalid mode. Should be direct_mapping, optimisation_and_conservation_1 or optimisation_and_conservation_2."
            )

    cleared_output_sequences = clear_output_sequences(output_sequences)

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

    output_list = []

    for name, input, output, identity_percentage in zip(
        map(lambda record: record.name, nucleotide_sequences),
        map(lambda record: str(record.seq), nucleotide_sequences),
        cleared_output_sequences,
        identity_percentages,
    ):
        # print(name, len(input) == len(output), len(output) - len(input))
        output_list.append(
            {
                "name": name,
                "input": input,
                "output": output,
                "identity_percentage": identity_percentage,
            }
        )

    return output_list
