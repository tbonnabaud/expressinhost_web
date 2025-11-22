import os
import random
import warnings
from dataclasses import dataclass
from typing import Iterable

import numpy as np
from Bio.PDB import DSSP, PDBParser, ShrakeRupley
from Bio.Seq import Seq

from ..logger import logger
from .codon_tables import ProcessedCodonTable, process_codon_table_from_file
from .exceptions import ExpressInHostError

AMINO_ACID_LETTER_MAPPING = {
    "G": "Gly",
    "P": "Pro",
    "A": "Ala",
    "V": "Val",
    "L": "Leu",
    "I": "Ile",
    "M": "Met",
    "C": "Cys",
    "F": "Phe",
    "Y": "Tyr",
    "W": "Trp",
    "H": "His",
    "K": "Lys",
    "R": "Arg",
    "Q": "Gln",
    "N": "Asn",
    "E": "Glu",
    "D": "Asp",
    "S": "Ser",
    "T": "Thr",
}


# Reference SASA values for each amino acid (in Å²) from Miller et al. 1987
# These are the SASA values for an amino acid in extended state of a protein (G-X-G)
REFERENCE_SASA = {
    "A": 113.0,
    "R": 241.0,
    "N": 158.0,
    "D": 151.0,
    "C": 140.0,
    "Q": 189.0,
    "E": 183.0,
    "G": 85.0,
    "H": 194.0,
    "I": 182.0,
    "L": 180.0,
    "K": 211.0,
    "M": 204.0,
    "F": 218.0,
    "P": 143.0,
    "S": 122.0,
    "T": 146.0,
    "W": 259.0,
    "Y": 229.0,
    "V": 160.0,
}


@dataclass
class Residue:
    residue_number: int
    amino_acid: str
    secondary_structure: str
    sasa: np.float64
    rsa: np.float64


@dataclass
class StructureInfos:
    name: str | None
    residue_list: list[Residue]


def extract_structure_infos(pdb_filename: str) -> StructureInfos | None:
    """
    Extract residue number, amino acid, secondary structure, SASA
    (solvent accessibility surface area) and RSA (relative SASA).
    """
    try:
        # Use a context manager to capture warnings
        with warnings.catch_warnings():
            # Treat UserWarnings as errors to handle DSSP parsing issues (e.g., malformed PDB)
            warnings.simplefilter("error", UserWarning)
            parser = PDBParser()
            structure = parser.get_structure("protein", pdb_filename)
            name = parser.get_header().get("name")
            model = structure[0]

            # DSSP for secondary structure
            dssp = DSSP(model, pdb_filename)

        # SASA values extraction
        sr = ShrakeRupley()
        sr.compute(structure, level="S")

        residue_list = []

        for key in dssp.keys():
            chain_id, res_id = key
            residue = model[chain_id][res_id]

            # Exclude heteroatoms and ensure it's an amino acid
            if res_id[0] == " " and "CA" in residue:
                resnum = res_id[1]
                aa = dssp[key][1]  # One-letter amino acid
                ss = dssp[key][2]  # Secondary structure code

                # computing RSA
                sasa = sum(atom.sasa for atom in residue if hasattr(atom, "sasa"))
                ref_sasa = REFERENCE_SASA.get(aa, 100.0)
                rsa = sasa / ref_sasa

                residue_list.append(
                    Residue(
                        residue_number=resnum,
                        amino_acid=aa,
                        secondary_structure=ss,
                        sasa=sasa,
                        rsa=rsa,
                    )
                )

        return StructureInfos(name=name, residue_list=residue_list)

    except Exception as exc:
        raise ExpressInHostError(str(exc))


def select_codon_from_table(
    amino_acid: str, host_codon_table: ProcessedCodonTable, rsa: float
) -> str:
    """
    Selects a codon for a given amino acid based on a host codon table and relative solvent accessibility (RSA).

    This function filters the codons corresponding to the specified amino acid from the host codon table
    and selects one codon based on a weighted probability that considers the RSA and the rank of each codon.

    Parameters:
        amino_acid (str): A single-letter code representing the amino acid.
        host_codon_table (ProcessedCodonTable): A processed codon table containing codon information.
        rsa (float): The relative solvent accessibility value.

    Returns:
        str: The selected codon as a string.

    Notes:
    - If the amino acid is Methionine ('M'), the function returns 'AUG' directly.
    - Else if the amino acid is Tryptophan ('W'), the function returns 'UGG' directly.
    - The probability of selecting a codon is calculated as:
      `Probability(codon) = RSA + Rank(codon) * (1 - 2 * RSA)`
    """
    # Methionine has only one possibility
    if amino_acid == "M":
        return "AUG"

    # Tryptophan has only one possibility
    elif amino_acid == "W":
        return "UGG"

    # RSA above 1 is potentially either at the beginning or at the end of the structure
    # In this case, set value to 1 to avoid negative weights
    if rsa > 1:
        rsa = 1

    # Keep rows with a specific amino-acid
    filtered_rows = [
        row
        for row in host_codon_table.values()
        # Use mapping to convert the single letter into to three letters
        if row.amino_acid == AMINO_ACID_LETTER_MAPPING[amino_acid]
    ]

    available_codons = [row.codon for row in filtered_rows]
    # Probability(codon) = RSA + Rank(codon) * (1 - 2 * RSA)
    codon_weights = [rsa + row.rank * (1 - 2 * rsa) for row in filtered_rows]

    return random.choices(available_codons, weights=codon_weights)[0]


def generate_mrna_from_residue_list(
    residue_list: Iterable[Residue],
    host_codon_table: ProcessedCodonTable,
) -> str:
    """
    Generates an mRNA sequence from a list of residues using a host codon table.

    This function iterates over a list of residues, selects a codon for each residue based on its
    amino acid and relative solvent accessibility (RSA), and concatenates the selected codons to form
    an mRNA sequence.

    Parameters:
        residue_list (Iterable[Residue]): An iterable of Residue objects, each representing an amino acid residue.
        host_codon_table (ProcessedCodonTable): A processed codon table containing codon information.

    Returns:
        str: The generated mRNA sequence as a string of concatenated codons.

    Notes:
    - The function uses the `select_codon_from_table` function to select a codon for each residue.
    - The RSA value of each residue is used to determine the probability of selecting each codon.
    """

    def generator():
        for residue in residue_list:
            yield select_codon_from_table(
                residue.amino_acid, host_codon_table, residue.rsa
            )

    return "".join(generator())


# Example usage
if __name__ == "__main__":
    pdb_file = "examples/AF-P37330-F1-model_v4.pdb"
    host_codon_table = process_codon_table_from_file(
        "codon_tables/Bacillus_subtilis.csv", 0.5
    )

    if os.path.exists(pdb_file):
        structure_infos = extract_structure_infos(pdb_file)

        if structure_infos is not None:
            input_aa = "".join(
                [residue.amino_acid for residue in structure_infos.residue_list]
            )
            mrna_sequence = generate_mrna_from_residue_list(
                structure_infos.residue_list, host_codon_table
            )
            output_aa = Seq(mrna_sequence).translate()
            assert input_aa == str(output_aa)

            print(f"Name: {structure_infos.name}\n")
            print(f"Amino-acid sequence:\n{input_aa}\n")
            print(f"mRNA sequence:\n{mrna_sequence}")

    else:
        print(f"PDB file not found: {pdb_file}")
