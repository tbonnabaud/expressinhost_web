def clear_output_sequences(output_sequences: list[str]) -> list[str]:
    return [seq.replace("-", "").replace("U", "T") for seq in output_sequences]


def compare_sequences(
    cleared_nucleotide_sequences: list[str], cleared_output_sequences: list[str]
) -> list:
    results = []

    for input_seq, output_seq in zip(
        cleared_nucleotide_sequences, cleared_output_sequences
    ):
        counter_similarity = 0

        for t in range(int(len(input_seq) / 3)):
            input_codon = input_seq[3 * t : 3 * t + 3].replace("U", "T")
            output_codon = output_seq[3 * t : 3 * t + 3]

            if input_codon == output_codon:
                counter_similarity += 1

        percentage = counter_similarity / (len(input_seq) / 3) * 100
        results.append(percentage)

    return results
