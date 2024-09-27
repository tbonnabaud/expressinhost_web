from pathlib import Path

import polars as pl

from src.codon_tables import process_codon_table_from_file
from src.constantes import CONSERVATION_THRESHOLD, SLOW_SPEED_THRESHOLD, WOBBLE_RATE
from src.schemas import TuningParameters
from src.utils import (
    find_organism_from_nucleotide_name,
    parse_sequences,
    read_text_file,
)


def get_output_sequences(content: str) -> list[str]:
    return content.split("\n")


def find_amino_acid_and_rank_and_speed(
    codon: str, table: list[dict]
) -> tuple[str, float] | tuple[None, None, None]:
    """Return the tuple corresponding to the given codon."""
    for row in table:
        if row["Codon"] == codon:
            return row["AA"], row["Rank"], row["Symbol_Speed"]

    return None, None, None


def generate_rank_lists(
    input_sequence_path: Path,
    output_sequence_path: Path,
    host_organism: str,
    results_path: Path,
):
    tuning_parameters = TuningParameters(
        wobble_rate=WOBBLE_RATE,
        slow_speed_threshold=SLOW_SPEED_THRESHOLD,
        conservation_threshold=CONSERVATION_THRESHOLD,
    )

    input_content = read_text_file(input_sequence_path)
    input_records = parse_sequences(input_content, "fasta")

    output_content = read_text_file(output_sequence_path)
    output_sequences = get_output_sequences(output_content)

    native_organism_list = [
        find_organism_from_nucleotide_name(record.name) for record in input_records
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

    host_codon_rows = host_codon_table.rows(named=True)

    for input_record, output_seq, native_codon_table in zip(
        input_records, output_sequences, native_codon_tables
    ):
        native_codon_rows = native_codon_table.rows(named=True)
        input_seq = input_record.seq.transcribe()
        output_seq = output_seq.replace("T", "U")

        print(input_record.id)
        print("-" * 10)

        input_rank_list = []
        input_symbol_speed_list = []
        output_rank_list = []

        for t in range(int(len(input_record.seq) / 3)):
            input_codon = input_seq[3 * t : 3 * t + 3]
            output_codon = output_seq[3 * t : 3 * t + 3]

            input_aa, input_rank, input_symbol_speed = (
                find_amino_acid_and_rank_and_speed(input_codon, native_codon_rows)
            )
            _, output_rank, _ = find_amino_acid_and_rank_and_speed(
                output_codon, host_codon_rows
            )

            if input_aa:
                input_rank_list.append(input_rank)
                output_rank_list.append(output_rank)
                input_symbol_speed_list.append(input_symbol_speed)

        print(len(input_rank_list), len(output_rank_list))
        df = pl.DataFrame(
            {
                "input_ranks": input_rank_list,
                "symbol_speed": input_symbol_speed_list,
                "output_ranks": output_rank_list,
            }
        )
        df.write_csv(results_path / f"{input_record.id}_ranks.csv")


def run(
    input_sequence_path: Path,
    output_sequence_path: Path,
    host_organism: str,
    results_path: Path,
):
    print(results_path)
    print("=" * 10)
    results_path.mkdir(parents=True, exist_ok=True)
    generate_rank_lists(
        input_sequence_path, output_sequence_path, host_organism, results_path
    )


if __name__ == "__main__":
    host_organism = "Escherichia_coli"

    run(
        Path("examples/Rad51_nucleotide.txt"),
        Path("output/Escherichia_coli/optimisation_and_conservation_1.txt"),
        host_organism,
        Path("results/new_software"),
    )

    run(
        Path("examples/Rad51_nucleotide.txt"),
        Path.home()
        / "workspace/HHU/expressinhost/Output/Example_Optim1_E_Coli/Optimisation_and_conservation_1.txt",
        host_organism,
        Path("results/old_software"),
    )

    # run(
    #     Path("examples/Rad51_nucleotide_bis.txt"),
    #     Path("output/Escherichia_coli/direct_mapping.txt"),
    #     host_organism,
    #     Path("results/new_software"),
    # )

    # run(
    #     Path("examples/Rad51_nucleotide_bis.txt"),
    #     Path.home() / "workspace/HHU/expressinhost/Direct_mapping.txt",
    #     host_organism,
    #     Path("results/old_software"),
    # )
