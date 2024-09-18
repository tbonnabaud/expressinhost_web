import argparse
import json
from functools import cache
from pathlib import Path

import polars as pl

from .codon_tables import process_raw_codon_table
from .constantes import TABLE_HEADERS
from .postprocessing import clear_output_sequences, compare_sequences
from .preprocessing import align_nucleotide_sequences, clear_nucleotide_sequences
from .simulations import (
    direct_mapping,
    optimisation_and_conservation_1,
    optimisation_and_conservation_2,
)
from .utils import (
    parse_alignments,
    parse_sequences,
    read_text_file,
    timeit,
    write_text_to_file,
)

TABLE_BASE_PATH = Path("codon_tables/")


def process_codon_table_from_file(organism_name: str) -> pl.DataFrame:
    table_df = pl.read_csv(
        TABLE_BASE_PATH / f"{organism_name}.txt",
        has_header=False,
        separator="\t",
        new_columns=TABLE_HEADERS,
    )
    processed_df = process_raw_codon_table(table_df)
    processed_df.write_csv(
        f"tmp/processed_tables/Processed_{organism_name}.csv",
        separator="\t",
    )

    return processed_df


@cache
def get_available_organism_list():
    return [file.name.replace(".txt", "") for file in TABLE_BASE_PATH.iterdir()]


def find_organism_from_nucleotide_name(name: str) -> str:
    for organism in get_available_organism_list():
        if organism.lower() in name.lower():
            return organism

    raise Exception(f"Organism not found for nucleotide: {name}.")


@timeit
def run_simulation(
    nucleotide_file_content: str,
    clustal_file_content: str | None,
    host_organism: str,
    mode: str,
):
    nucleotide_sequences = parse_sequences(nucleotide_file_content, "fasta")
    cleared_nucleotide_sequences = clear_nucleotide_sequences(nucleotide_sequences)

    native_organism_list = [
        find_organism_from_nucleotide_name(record.name)
        for record in nucleotide_sequences
    ]

    native_codon_tables = [
        process_codon_table_from_file(name) for name in native_organism_list
    ]

    host_codon_table = process_codon_table_from_file(host_organism)

    if mode == "direct_mapping":
        output_sequences = direct_mapping(
            cleared_nucleotide_sequences, native_codon_tables, host_codon_table
        )

    else:
        if clustal_file_content is None:
            raise Exception("Clustal file is required.")

        clustal_sequences = parse_sequences(clustal_file_content, "clustal")

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

    with open(output_path / f"{mode}_identity_percentages.json", "w") as f:
        name_result_mapping = dict(zip(nucleotide_names, identity_percentages))
        json.dump(name_result_mapping, f, indent=4)


def main(
    nucleotide_file_path: Path,
    clustal_file_path: Path | None,
    host_organism: str,
    mode: str,
):
    tmp_dirpath = Path("tmp/")
    (tmp_dirpath / "processed_tables").mkdir(parents=True, exist_ok=True)

    # Remove existing files
    for path in tmp_dirpath.iterdir():
        if path.is_file():
            path.unlink()

    nucleotide_file_content = read_text_file(nucleotide_file_path)

    if clustal_file_path is None:
        run_simulation(nucleotide_file_content, None, host_organism, mode)

    else:
        clustal_file_content = read_text_file(clustal_file_path)
        run_simulation(
            nucleotide_file_content, clustal_file_content, host_organism, mode
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ExpressInHost CLI")

    parser.add_argument(
        "--nucleotide-file",
        type=str,
        required=True,
        help="Path to nucleotide file (FASTA format).",
    )

    parser.add_argument(
        "--clustal-file",
        type=str,
        default=None,
        help='Path to clustal file. Header should begin by "CLUSTAL".',
    )

    parser.add_argument(
        "--host",
        type=str,
        required=True,
        choices=[
            "Arabidopsis_thaliana",
            "Bacillus_subtilis",
            "Caenorhabditis_elegans",
            "Danio_rerio",
            "Drosophila_melanogaster",
            "Escherichia_coli",
            "Gallus_gallus",
            "Homo_sapiens",
            "Komagataella_pastoris",
            "Methanocaldococcus_jannaschii",
            "Saccharomyces_cerevisiae",
            "Staphylococcus_aureus",
            "Xenopus_laevis",
        ],
        help="Host organism.",
    )

    parser.add_argument(
        "--mode",
        type=str,
        required=True,
        choices=[
            "direct_mapping",
            "optimisation_and_conservation_1",
            "optimisation_and_conservation_2",
        ],
        help="Simulation mode.",
    )

    args = parser.parse_args()

    main(args.nucleotide_file, args.clustal_file, args.host, args.mode)
