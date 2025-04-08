def clear_output_sequence(sequence: str) -> str:
    """Remove dashes."""
    return sequence.replace("-", "")


def compute_similarity(input_sequence: str, output_sequence: str) -> float:
    """
    Compute the similarity between input and output sequences by comparing codons.
    """
    counter_similarity = 0
    codon_nb = int(len(input_sequence) / 3)

    for t in range(codon_nb):
        input_codon = input_sequence[3 * t : 3 * t + 3].replace("U", "T")
        output_codon = output_sequence[3 * t : 3 * t + 3]

        if input_codon == output_codon:
            counter_similarity += 1

    return counter_similarity / codon_nb * 100
