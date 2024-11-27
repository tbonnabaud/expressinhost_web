import csv
from dataclasses import dataclass

import numpy as np

from .constantes import AMINO_ACID_LIST, TABLE_BASE_PATH


@dataclass(slots=True)
class ProcessedCodonTableRow:
    amino_acid: str
    anticodon: str
    codon: str
    trna_gcn: float
    corresp_codon: str
    wobble_rate: float
    rank: float
    speed: float
    symbol_speed: str


@dataclass(slots=True)
class ProcessedCodonTable:
    indexed_rows: dict[str, ProcessedCodonTableRow]

    def __getitem__(self, key: str):
        return self.indexed_rows[key]

    def get(self, key: str):
        return self.indexed_rows.get(key)

    def values(self):
        return self.indexed_rows.values()


def process_raw_codon_table(
    raw_codon_table: list[dict], slow_speed_threshold: float
) -> ProcessedCodonTable:
    """Any raw data table contains:
    - column 1: amino acid (amino_acid)
    - column 2: anti-codon (anticodon)
    - column 3: codon
    - column 4: tRNA GCN (trna_gcn)
    - column 5: codon used in case of wobble (corresp_codon)
    - column 6: wobble rate (wobble_rate)

    All input raw tables contain 61 lines (for 61 codons).

    Only column 4 and 5 need to be adapted to your specific organism.
    To find the tRNA GCN for your organism of interest,
    use databases such as "GtRNAdb: Genomic tRNA Database".

    Args:
        raw_codon_table (list[dict]): raw table

    Returns:
        ProcessedCodonTable: indexed processed table with three more columns (rank, speed and symbol_speed)
    """
    # Column 1 of the input raw table
    amino_acid_col = [row["amino_acid"] for row in raw_codon_table]
    # Column 2 of the input raw table
    anticodon_col = [row["anticodon"] for row in raw_codon_table]
    # Column 3 of the input raw table
    codon_col = [row["codon"] for row in raw_codon_table]
    # Column 4 of the input raw table
    gcn_col = np.array([float(row["trna_gcn"]) for row in raw_codon_table])
    # Column 5 of the input raw table
    corresp_codon_col = [row["corresp_codon"] for row in raw_codon_table]
    # Column 6 of the input raw table
    wobble_rate_col = np.array([float(row["wobble_rate"]) for row in raw_codon_table])

    # Colmun 7 to be added to the output processed table
    rank_col = np.zeros(61, dtype=np.float64)
    # Column 8 to be added to the output processed table
    speed_col = np.zeros(61, dtype=np.float64)
    # Column 9 to be added to the output processed table
    symbol_speed_col = ["0" for _ in range(61)]

    # This index allows to know with which amino acid of the list we are dealing (the first is "Ala" the second is "Arg"...)
    aa_num = 0
    # For each amino acid, it counts the number of codons that have been processed
    cpt_codons_aa = 0
    # Highest GCN for each amino acid
    top_gcn_aa = 0.0
    # Lowest GCN for each amino acid
    last_gcn_aa = 1000000.0
    # Stores the sum of all GCN of the output processed table
    gcn_tot = 0.0

    # We complement the input raw table for the codons with GCN = 0 and we count the total number of GCN, including the complemented ones
    k = 0
    while k < 61:
        if amino_acid_col[k] == AMINO_ACID_LIST[aa_num]:
            # If no GCN for that codon
            if gcn_col[k] == 0:
                # In all codons of the input raw table
                for t in range(61):
                    # We search its corresponding codon
                    if codon_col[t] == corresp_codon_col[k]:
                        # And calculate its GCN using the GCN of the codon it wobbles with, and the wobbling rate
                        gcn_col[k] = gcn_col[t] - wobble_rate_col[k] * gcn_col[t]
                        break

            # If the GCN is higher than the top one
            if gcn_col[k] > top_gcn_aa:
                top_gcn_aa = gcn_col[k]

            # If the GCN is lower than the lowest one
            if gcn_col[k] < last_gcn_aa:
                last_gcn_aa = gcn_col[k]

            # One more codon for that amino acid has been processed
            cpt_codons_aa += 1
            gcn_tot = gcn_tot + gcn_col[k]

        elif amino_acid_col[k] == AMINO_ACID_LIST[aa_num + 1]:
            # If more than one codon code for that amino acid
            if cpt_codons_aa > 1:
                # For all codons of the amino acid we are currently dealing with
                for t in range(cpt_codons_aa):
                    if top_gcn_aa != last_gcn_aa:
                        # The "RANK" of the codon is calculated using its GCN, the highest and the lowest GCN of that amino acid
                        rank_col[k - cpt_codons_aa + t] = (
                            gcn_col[k - cpt_codons_aa + t] - last_gcn_aa
                        ) / (top_gcn_aa - last_gcn_aa)
                    else:
                        # If the top and lowest codons have the same GCN, the "RANK" of the codon is 1
                        rank_col[k - cpt_codons_aa + t] = 1.0

            # If only one codon codes for that amino acid
            if cpt_codons_aa == 1:
                rank_col[k - cpt_codons_aa] = 1

            # Switch to the next amino acid of the amino acid list
            aa_num += 1
            # Reset the counter of codons for each amino acids
            cpt_codons_aa = 0
            # Reset the value for highest GCN
            top_gcn_aa = 0.0
            # Reset value for lowest GCN
            last_gcn_aa = 1000000.0
            # Step one line back such that this line itself will be processed by entering the first "If condition" of the "For loop"
            k -= 1

        if k == 60:
            # If more than one codon code for that amino acid
            if cpt_codons_aa > 1:
                # For all codons of the amino acid we are currently dealing with
                for t in range(cpt_codons_aa):
                    if top_gcn_aa != last_gcn_aa:
                        # The "RANK" of the codon is calculated using its GCN, the highest and the lowest GCN of that amino acid
                        rank_col[k + 1 - cpt_codons_aa + t] = (
                            gcn_col[k + 1 - cpt_codons_aa + t] - last_gcn_aa
                        ) / (top_gcn_aa - last_gcn_aa)
                    else:
                        # If the top and lowest codons have the same GCN, the "RANK" of the codon is 1
                        rank_col[k + 1 - cpt_codons_aa + t] = 1.0

            # If only one codon codes for that amino acid
            if cpt_codons_aa == 1:
                rank_col[k + 1 - cpt_codons_aa] = 1

        k += 1

    # Calculate the "SPEED" and the average "SPEED" over the entire table
    speed_col = gcn_col / gcn_tot
    average_speed = speed_col.mean()

    # Search the lowest "SPEED" over the entire table
    lowest_speed = speed_col.min()

    # Tag the codons of low "SPEED" over the entire table
    for k in range(61):
        # If codon's SPEED is below the threshold set by the slow_speed_threshold.
        # That threshold is a limit SPEED value, independent from the number of codons that can fall in that category.
        if speed_col[k] < lowest_speed + (
            slow_speed_threshold * (average_speed - lowest_speed)
        ):
            symbol_speed_col[k] = "S"

    return ProcessedCodonTable(
        indexed_rows={
            row[2]: ProcessedCodonTableRow(*row)
            for row in zip(
                amino_acid_col,
                anticodon_col,
                codon_col,
                gcn_col,
                corresp_codon_col,
                wobble_rate_col,
                rank_col,
                speed_col,
                symbol_speed_col,
            )
        }
    )


def process_codon_table_from_file(
    codon_table_name: str, slow_speed_threshold: float
) -> ProcessedCodonTable:
    with open(TABLE_BASE_PATH / f"{codon_table_name}.csv") as file:
        reader = csv.DictReader(file, delimiter="\t")

        return process_raw_codon_table(list(reader), slow_speed_threshold)
