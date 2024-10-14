import argparse
import re
from pathlib import Path

from .core.codon_tables import process_codon_table_from_file
from .core.constantes import CONSERVATION_THRESHOLD, SLOW_SPEED_THRESHOLD
from .core.sequence_tuning import run_tuning
from .core.utils import find_organism_from_nucleotide_name, read_text_file
from .schemas import TuningParameters


def create_and_clean_dirs():
    tmp_dirpath = Path("tmp/")
    (tmp_dirpath / "processed_tables").mkdir(parents=True, exist_ok=True)
    # Create output directory
    Path("output/").mkdir(exist_ok=True)

    # Remove existing files
    for path in tmp_dirpath.iterdir():
        if path.is_file():
            path.unlink()


def run_as_cli(
    nucleotide_file_path: Path,
    clustal_file_path: Path | None,
    host_codon_table_name: str,
    mode: str,
):
    nucleotide_file_content = read_text_file(nucleotide_file_path)

    tuning_parameters = TuningParameters(
        slow_speed_threshold=SLOW_SPEED_THRESHOLD,
        conservation_threshold=CONSERVATION_THRESHOLD,
    )

    sequence_names = re.findall(
        r"^\> ?([\w ]*\w)", nucleotide_file_content, re.MULTILINE
    )

    native_codon_table_names = [
        find_organism_from_nucleotide_name(name) for name in sequence_names
    ]

    native_codon_tables = [
        process_codon_table_from_file(name, tuning_parameters.slow_speed_threshold)
        for name in native_codon_table_names
    ]

    host_codon_table = process_codon_table_from_file(
        host_codon_table_name,
        tuning_parameters.slow_speed_threshold,
    )

    if clustal_file_path is None:
        run_tuning(
            nucleotide_file_content,
            None,
            native_codon_tables,
            host_codon_table,
            mode,
            tuning_parameters,
        )

    else:
        clustal_file_content = read_text_file(clustal_file_path)
        run_tuning(
            nucleotide_file_content,
            clustal_file_content,
            native_codon_tables,
            host_codon_table,
            mode,
            tuning_parameters,
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

    run_as_cli(args.nucleotide_file, args.clustal_file, args.host, args.mode)
