from typing import Iterable


def clear_output_sequences(output_sequences: Iterable[str]):
    for seq in output_sequences:
        yield seq.replace("-", "").replace("U", "T")


def compare_sequences(
    cleared_nucleotide_sequences: Iterable[str], cleared_output_sequences: Iterable[str]
):
    for input_seq, output_seq in zip(
        cleared_nucleotide_sequences, cleared_output_sequences
    ):
        counter_similarity = 0

        for t in range(int(len(input_seq) / 3)):
            input_codon = input_seq[3 * t : 3 * t + 3].replace("U", "T")
            output_codon = output_seq[3 * t : 3 * t + 3]

            if input_codon == output_codon:
                counter_similarity += 1

        yield counter_similarity / (len(input_seq) / 3) * 100
