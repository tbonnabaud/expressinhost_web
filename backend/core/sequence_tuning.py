import random

from .checks import check_amino_acido_conservation, check_nucleotides_clustal_identity
from .codon_tables import ProcessedCodonTable
from .exceptions import NoAminoAcidConservation, NoIdenticalSequencesError
from .postprocessing import (
    clear_output_sequence,
    compute_similarity,
)
from .preprocessing import align_nucleotide_sequences, dna_to_rna_sequences
from .utils import (  # write_text_to_file,
    get_clustal_symbol_sequence,
    parse_sequences,
    timeit,
)


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
                    # Add stop codon
                    new_line += input_codon

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
                    # Add stop codon
                    new_line += input_codon

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
                    # Add stop codon
                    new_line += input_codon

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


def get_sequence_profiles(sequence: str, codon_table: ProcessedCodonTable):
    """Get speed and rank profiles."""
    speed_profile = []
    rank_profile = []

    for t in range(int(len(sequence) / 3)):
        codon = sequence[3 * t : 3 * t + 3]
        row = codon_table.get(codon)

        if row:
            speed_profile.append(row.speed)
            rank_profile.append(row.rank)

    return {"speed": speed_profile, "rank": rank_profile}


class SequenceTuner:
    def __init__(
        self,
        nucleotide_file_content: str,
        clustal_file_content: str | None,
        native_codon_tables: list[ProcessedCodonTable],
        host_codon_table: ProcessedCodonTable,
        # mode: str,
        # conservation_threshold: float | None,
    ):
        self.nucleotide_records = parse_sequences(nucleotide_file_content, "fasta")
        self.clustal_records = (
            parse_sequences(clustal_file_content, "clustal")
            if clustal_file_content
            else None
        )
        self.symbol_sequence = (
            get_clustal_symbol_sequence(clustal_file_content)
            if clustal_file_content
            else None
        )
        self.native_codon_tables = native_codon_tables
        self.host_codon_table = host_codon_table

        if clustal_file_content:
            self.ensure_sequence_matching()

    def ensure_sequence_matching(self):
        """Ensure sequences are the same in the two files."""
        for nucleotide_record, clustal_record in zip(
            self.nucleotide_records, self.clustal_records
        ):
            if not check_nucleotides_clustal_identity(
                nucleotide_record, clustal_record
            ):
                raise NoIdenticalSequencesError(
                    f"Sequences {nucleotide_record.name} (FASTA) and {clustal_record.name} (CLUSTAL) are not identical."
                    "Check their value in the two files and check their order."
                )

    def process(self, mode: str, conservation_threshold: float | None) -> list[str]:
        cleared_nucleotide_sequences = dna_to_rna_sequences(self.nucleotide_records)

        if mode == "direct_mapping":
            return direct_mapping(
                cleared_nucleotide_sequences,
                self.native_codon_tables,
                self.host_codon_table,
            )

        else:
            if self.nucleotide_records is None:
                raise Exception("Clustal file is required.")

            aligned_nucleotide_sequences = align_nucleotide_sequences(
                self.clustal_records, cleared_nucleotide_sequences
            )

            if mode == "optimisation_and_conservation_1":
                return optimisation_and_conservation_1(
                    aligned_nucleotide_sequences,
                    self.symbol_sequence,
                    self.native_codon_tables,
                    self.host_codon_table,
                )

            elif mode == "optimisation_and_conservation_2":
                return optimisation_and_conservation_2(
                    aligned_nucleotide_sequences,
                    self.symbol_sequence,
                    self.native_codon_tables,
                    self.host_codon_table,
                    conservation_threshold,
                )

            else:
                raise Exception(
                    "Invalid mode. Should be direct_mapping, optimisation_and_conservation_1 or optimisation_and_conservation_2."
                )

    def postprocess(self, output_sequences: list[str]) -> list[dict]:
        output_list = []

        for input_record, output_sequence, native_codon_table in zip(
            self.nucleotide_records,
            output_sequences,
            self.native_codon_tables,
        ):
            cleared_output_sequence = clear_output_sequence(output_sequence)

            # Ensure input and ouput nucleotide sequence have the same amino-acid sequence
            if not check_amino_acido_conservation(
                input_record, cleared_output_sequence
            ):
                raise NoAminoAcidConservation(
                    "Amino acid sequences from input and output are supposed to be the same."
                )

            input_sequence = str(input_record.seq)
            identity_percentage = compute_similarity(
                input_sequence, cleared_output_sequence
            )
            output_list.append(
                {
                    "name": input_record.name,
                    "input": input_sequence,
                    "output": cleared_output_sequence,
                    "identity_percentage": identity_percentage,
                    "input_profiles": get_sequence_profiles(
                        input_sequence, native_codon_table
                    ),
                    "output_profiles": get_sequence_profiles(
                        cleared_output_sequence, self.host_codon_table
                    ),
                }
            )

        return output_list


@timeit
def tune_sequences(
    nucleotide_file_content: str,
    clustal_file_content: str | None,
    native_codon_tables: list[ProcessedCodonTable],
    host_codon_table: ProcessedCodonTable,
    mode: str,
    conservation_threshold: float | None,
) -> list[dict]:
    # Preprocessing
    nucleotide_records = parse_sequences(nucleotide_file_content, "fasta")
    cleared_nucleotide_sequences = dna_to_rna_sequences(nucleotide_records)

    # Processing
    if mode == "direct_mapping":
        output_sequences = direct_mapping(
            cleared_nucleotide_sequences, native_codon_tables, host_codon_table
        )

    else:
        if clustal_file_content is None:
            raise Exception("Clustal file is required.")

        clustal_records = parse_sequences(clustal_file_content, "clustal")

        # Ensure sequences are the same in the two files
        for nucleotide_record, clustal_record in zip(
            nucleotide_records, clustal_records
        ):
            if not check_nucleotides_clustal_identity(
                nucleotide_record, clustal_record
            ):
                raise NoIdenticalSequencesError(
                    f"Sequences {nucleotide_record.name} (FASTA) and {clustal_record.name} (CLUSTAL) are not identical."
                    "Check their value in the two files and check their order."
                )

        # write_text_to_file(
        #     "\n".join([str(r.seq) for r in clustal_sequences]),
        #     "tmp/modif_sequences_2.txt",
        # )

        aligned_nucleotide_sequences = align_nucleotide_sequences(
            clustal_records, cleared_nucleotide_sequences
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

    output_list = []

    # Postprocessing
    for input_record, output_sequence, native_codon_table in zip(
        nucleotide_records,
        output_sequences,
        native_codon_tables,
    ):
        cleared_output_sequence = clear_output_sequence(output_sequence)

        # Ensure input and ouput nucleotide sequence have the same amino-acid sequence
        if not check_amino_acido_conservation(input_record, cleared_output_sequence):
            raise NoAminoAcidConservation(
                "Amino acid sequences from input and output are supposed to be the same."
            )

        input_sequence = str(input_record.seq)
        identity_percentage = compute_similarity(
            input_sequence, cleared_output_sequence
        )
        output_list.append(
            {
                "name": input_record.name,
                "input": input_sequence,
                "output": cleared_output_sequence,
                "identity_percentage": identity_percentage,
                "input_profiles": get_sequence_profiles(
                    input_sequence, native_codon_table
                ),
                "output_profiles": get_sequence_profiles(
                    cleared_output_sequence, host_codon_table
                ),
            }
        )

    return output_list
