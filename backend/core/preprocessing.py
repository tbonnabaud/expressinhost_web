import re

from Bio.SeqRecord import SeqRecord


def clear_nucleotide_sequences(nucleotide_sequences: list[SeqRecord]) -> list[str]:
    stop_codon_pattern = re.compile(r"(UAA|UAG|UGA)$", re.IGNORECASE)
    cleared_nucleotide_sequences = []

    for record in nucleotide_sequences:
        # Transcribe DNA to RNA (replace T by U)
        seq = record.seq.transcribe()
        # Remove UAA|UAG|UGA stop codons
        str_seq = stop_codon_pattern.sub("", str(seq))

        cleared_nucleotide_sequences.append(str_seq)

    # write_text_to_file(
    #     "\n".join(cleared_nucleotide_sequences), "tmp/modif_sequences_4.txt"
    # )

    return cleared_nucleotide_sequences


def align_nucleotide_sequences(
    clustal_sequences: list[SeqRecord], cleared_nucleotide_sequences: list[str]
) -> list[str]:
    aligned_nucleotide_sequences = []

    for clustal_record, nucleo_seq in zip(
        clustal_sequences, cleared_nucleotide_sequences
    ):
        new_line = ""

        m = 0

        for t in range(len(clustal_record.seq)):
            if clustal_record.seq[t] == "-":
                new_line += "---"

            else:
                new_line += nucleo_seq[3 * m : 3 * m + 3]
                m += 1

        aligned_nucleotide_sequences.append(new_line + "UAA")

    # write_text_to_file(
    #     "\n".join(aligned_nucleotide_sequences), "tmp/modif_sequences_5.txt"
    # )

    return aligned_nucleotide_sequences
