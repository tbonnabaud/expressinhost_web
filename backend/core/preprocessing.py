from Bio.SeqRecord import SeqRecord


def dna_to_rna_sequences(nucleotide_sequences: list[SeqRecord]) -> list[str]:
    return [str(record.seq.transcribe()) for record in nucleotide_sequences]


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

        # Keep the pre-tuned stop codon
        aligned_nucleotide_sequences.append(new_line + nucleo_seq[-3:])

    # write_text_to_file(
    #     "\n".join(aligned_nucleotide_sequences), "tmp/modif_sequences_5.txt"
    # )

    return aligned_nucleotide_sequences
