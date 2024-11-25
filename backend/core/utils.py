import time
from functools import cache
from io import StringIO
from pathlib import Path
from typing import Literal

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

from ..logger import logger
from .constantes import TABLE_BASE_PATH
from .exceptions import ClustalFormatError, FastaFormatError


def timeit(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        elapsed_time = end_time - start_time
        logger.info(f"Elapsed time for '{func.__name__}': {elapsed_time:.4f} seconds")

        return result  # Return the result of the function

    return wrapper


def read_text_file(path: Path | str) -> str:
    with open(path) as f:
        return f.read()


def write_text_to_file(content: str, path: Path | str):
    with open(path, "w") as f:
        f.write(content)


def parse_sequences(
    raw_content: str, format: Literal["fasta", "clustal"]
) -> list[SeqRecord]:
    try:
        return [record for record in SeqIO.parse(StringIO(raw_content), format)]

    except Exception:
        if format == "clustal":
            raise ClustalFormatError(
                "Fail to parse file. Ensure header start with CLUSTAL and file is correctly formatted."
            )

        else:
            raise FastaFormatError(
                "Fail to parse file. Ensure file is correctly formatted."
            )


def parse_alignments(
    raw_content: str, format: Literal["fasta", "clustal"]
) -> list[MultipleSeqAlignment]:
    try:
        return [record for record in AlignIO.parse(StringIO(raw_content), format)]

    except Exception:
        if format == "clustal":
            ClustalFormatError(
                "Fail to parse file. Ensure header start with CLUSTAL and file is correctly formatted."
            )

        else:
            FastaFormatError("Fail to parse file. Ensure file is correctly formatted.")


def get_clustal_symbol_sequence(clustal_file_content: str) -> str:
    clustal_alignments = parse_alignments(clustal_file_content, "clustal")[0]
    return clustal_alignments.column_annotations.get("clustal_consensus")


@cache
def get_available_organism_list():
    return [file.name.replace(".csv", "") for file in TABLE_BASE_PATH.iterdir()]


def find_organism_from_nucleotide_name(name: str) -> str:
    """Guess organism name from nucleotide name."""
    for organism in get_available_organism_list():
        if organism.lower() in name.lower():
            return organism

    raise Exception(f"Organism not found for nucleotide: {name}.")
