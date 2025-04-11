import os
import random

import pandas as pd
from Bio.PDB import DSSP, PDBParser, ShrakeRupley

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


def extract_structure_info(pdb_filename: str):
    """
    Extract residue number, amino acid, secondary structure, SASA
    (solvent accessibility surface area) and RSA (relative SASA).
    """
    try:
        parser = PDBParser()
        structure = parser.get_structure("protein", pdb_filename)
        model = structure[0]

        # DSSP for secondary structure
        dssp = DSSP(model, pdb_filename)

        # SASA values extraction
        sr = ShrakeRupley()
        sr.compute(structure, level="S")

        # storage lists
        residue_numbers = []
        amino_acids = []
        secondary_structures = []
        sasa_values = []
        rsa_values = []

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
                rsa = (sasa / ref_sasa) * 100

                # Append to lists
                residue_numbers.append(resnum)
                amino_acids.append(aa)
                secondary_structures.append(ss)
                sasa_values.append(sasa)
                rsa_values.append(rsa)

        # Create DataFrame
        return pd.DataFrame(
            {
                "Residue Number": residue_numbers,
                "Amino Acid": amino_acids,
                "Secondary Structure": secondary_structures,
                "SASA": sasa_values,
                "RSA": rsa_values,
            }
        )

    except Exception as e:
        print(f"Error processing {pdb_filename}: {str(e)}")
        return None


def select_codon_from_table(
    amino_acid: str, host_codon_table: ProcessedCodonTable, is_slow: bool
):
    # Keep rows with a specific amino-acid
    filtered_rows = [
        row
        for row in host_codon_table.values()
        if row.amino_acid == AMINO_ACID_LETTER_MAPPING[amino_acid]
    ]
    # total_gcn = sum(map(lambda x: x.trna_gcn, filtered_rows))

    slow_rows = [row for row in filtered_rows if row.rank < 0.5]
    fast_rows = [row for row in filtered_rows if row.rank >= 0.5]

    if is_slow:
        if slow_rows:
            available_codons = [row.codon for row in slow_rows]
            codon_weights = [(1 - row.rank) ** 2 for row in slow_rows]

        # Else return the slowest of the fast codons
        else:
            return min(fast_rows, key=lambda row: row.rank).codon

    else:
        if fast_rows:
            available_codons = [row.codon for row in fast_rows]
            codon_weights = [row.rank**2 for row in fast_rows]

        # Else return the fastest of the slow codons
        else:
            return max(slow_rows, key=lambda row: row.rank).codon

    return random.choices(available_codons, weights=codon_weights)[0]


# Example usage
if __name__ == "__main__":
    pdb_file = "examples/AF-P37330-F1-model_v4.pdb"
    host_codon_table = process_codon_table_from_file(
        "codon_tables/Bacillus_subtilis.csv", 0.5
    )

    if os.path.exists(pdb_file):
        df = extract_structure_info(pdb_file)

        if df is not None:
            # print(df.describe())

            for amino_acid, structure, rsa in zip(
                df["Amino Acid"], df["Secondary Structure"], df["RSA"]
            ):
                # print(amino_acid, structure, rsa)
                # RSA > to 25 means outside, so considered as slow
                is_slow = rsa > 25
                codon = select_codon_from_table(amino_acid, host_codon_table, is_slow)
                print(f"{codon}, {is_slow=}")

            # print(df)
    else:
        print(f"PDB file not found: {pdb_file}")
