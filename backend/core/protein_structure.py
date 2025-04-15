import os
import random
from dataclasses import dataclass
from typing import Iterable

import numpy as np
from Bio.PDB import DSSP, PDBParser, ShrakeRupley
from Bio.Seq import Seq

from .codon_tables import ProcessedCodonTable, process_codon_table_from_file

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

    except Exception as e:
        print(f"Error processing {pdb_filename}: {str(e)}")
        return None


def select_codon_from_table(
    amino_acid: str, host_codon_table: ProcessedCodonTable, is_slow: bool
) -> str:
    """
    Selects a codon for a given amino acid from a codon table based on the desired speed.

    This function filters the codon table to find codons that correspond to the specified
    amino acid. It then selects a codon based on whether a slow or fast codon is desired.
    If no codons are available for the desired speed, it selects the slowest fast codon
    or the fastest slow codon as a fallback.

    Parameters:
        amino_acid (str): The single-letter code of the amino acid for which to select a codon.
        host_codon_table (ProcessedCodonTable): A table containing codon information, including
                                                amino acid, codon, and rank.
        is_slow (bool): A flag indicating whether to select a slow codon (True) or a fast codon (False).

    Returns:
        str: The selected codon for the specified amino acid.

    Notes:
    - The function uses a rank to determine the speed of a codon, where a lower rank indicates
      a slower codon and a higher rank indicates a faster codon.
    - The function uses weighted random selection to choose a codon based on its rank to not exhaust
      the most abundant tRNA.
    - If no codons are available for the desired speed, the function falls back to the slowest
      fast codon or the fastest slow codon.
    """
    # Keep rows with a specific amino-acid
    filtered_rows = [
        row
        for row in host_codon_table.values()
        # Use mapping to convert the single letter into to three letters
        if row.amino_acid == AMINO_ACID_LETTER_MAPPING[amino_acid]
    ]

    slow_rows = [row for row in filtered_rows if row.rank < 0.5]
    fast_rows = [row for row in filtered_rows if row.rank >= 0.5]

    if is_slow:
        if slow_rows:
            available_slow_codons = [row.codon for row in slow_rows]
            slow_codon_weights = [(1 - row.rank) ** 2 for row in slow_rows]

            return random.choices(available_slow_codons, weights=slow_codon_weights)[0]

        # Else return the slowest of the fast codons
        else:
            return min(fast_rows, key=lambda row: row.rank).codon

    else:
        if fast_rows:
            available_fast_codons = [row.codon for row in fast_rows]
            fast_codon_weights = [row.rank**2 for row in fast_rows]

            return random.choices(available_fast_codons, weights=fast_codon_weights)[0]

        # Else return the fastest of the slow codons
        else:
            return max(slow_rows, key=lambda row: row.rank).codon


def generate_mrna_from_residue_list(
    residue_list: Iterable[Residue],
    host_codon_table: ProcessedCodonTable,
    rsa_threshold: float = 0.25,
) -> str:
    def generator():
        for residue in residue_list:
            # RSA > rsa_threshold means outside, so considered as slow
            is_slow = residue.rsa > rsa_threshold

            yield select_codon_from_table(residue.amino_acid, host_codon_table, is_slow)

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
                structure_infos.residue_list, host_codon_table, 0.25
            )
            output_aa = Seq(mrna_sequence).translate()
            assert input_aa == str(output_aa)

            print(f"Name: {structure_infos.name}\n")
            print(f"Amino-acid sequence:\n{input_aa}\n")
            print(f"mRNA sequence:\n{mrna_sequence}")

    else:
        print(f"PDB file not found: {pdb_file}")
