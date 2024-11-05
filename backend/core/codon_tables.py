from dataclasses import dataclass

import polars as pl

from .constantes import AMINO_ACID_LIST, TABLE_BASE_PATH


@dataclass
class ProcessedCodonTableRow:
    amino_acid: str
    anticodon: str
    codon: str
    trna_gcn: str
    corresp_codon: str
    wobble_rate: float
    rank: float
    speed: float
    symbol_speed: str


def process_raw_codon_table(
    raw_codon_table: pl.DataFrame, slow_speed_threshold: float
) -> list[dict]:
    """Any raw data table contains:
    - column 1: amino acid (AA)
    - column 2: anti-codon (tRNA)
    - column 3: codon (Codon)
    - column 4: tRNA GCN (GCN)
    - column 5: codon used in case of wobble (Corresp_codon)
    - column 6: wobble rate

    All input raw tables contain 61 lines (for 61 codons).

    Only column 4 and 5 need to be adapted to your specific organism.
    To find the tRNA GCN for your organism of interest,
    use databases such as "GtRNAdb: Genomic tRNA Database".

    Args:
        raw_codon_table (pl.DataFrame): raw table

    Returns:
        pl.DataFrame: processed table with three more columns (Rank, Speed and Symbol_Speed)
    """
    # Column 1 of the input raw table
    amino_acid_col = raw_codon_table["amino_acid"]
    # Column 2 of the input raw table
    anticodon_col = raw_codon_table["anticodon"]
    # Column 3 of the input raw table
    codon_col = raw_codon_table["codon"]
    # Column 4 of the input raw table
    gcn_col = raw_codon_table["trna_gcn"].cast(pl.Float64)
    # Column 5 of the input raw table
    corresp_codon_col = raw_codon_table["corresp_codon"]
    # Column 6 of the input raw table
    wobble_rate_col = raw_codon_table["wobble_rate"]

    # Colmun 7 to be added to the output processed table
    rank_col = pl.zeros(61, dtype=pl.Float64, eager=True).alias("rank")
    # Column 8 to be added to the output processed table
    speed_col = pl.zeros(61, dtype=pl.Float64, eager=True).alias("speed")
    # Column 9 to be added to the output processed table
    symbol_speed_col = pl.repeat("0", 61, dtype=pl.String, eager=True).alias(
        "symbol_speed"
    )

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
    speed_col = gcn_col.clone().alias("speed") / gcn_tot
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

    df = pl.DataFrame(
        [
            amino_acid_col,
            anticodon_col,
            codon_col,
            gcn_col,
            corresp_codon_col,
            wobble_rate_col,
            rank_col,
            speed_col,
            symbol_speed_col,
        ]
    )

    return df.rows(named=True)


def process_codon_table_from_file(
    codon_table_name: str, slow_speed_threshold: float
) -> pl.DataFrame:
    table_df = pl.read_csv(
        TABLE_BASE_PATH / f"{codon_table_name}.csv",
        has_header=True,
        separator="\t",
    )
    processed_df = process_raw_codon_table(table_df, slow_speed_threshold)
    # processed_df.write_csv(
    #     f"tmp/processed_tables/Processed_{codon_table_name}.csv",
    #     separator="\t",
    # )

    return processed_df
