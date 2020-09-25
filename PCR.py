from time import time
from typing import List, Tuple

from matplotlib import pyplot as plt

from GenBank import Sequence
from Utils import (clean_empty_sequences, clean_empty_tuples,
                   distance_between_primers, generate_fall_off_rate)


def PCR(segment: Tuple[Sequence, Sequence], primers: Tuple[Sequence, Sequence], num_cycles: int = 20, fall_off_noise: int = None) -> List[Tuple[Sequence, Sequence]]:
    """Simulate Polymerase Chain Reaction

    Args:
        segment (Tuple[Sequence, Sequence]): Segment to copy
        primers (Tuple[Sequence, Sequence]): Primers on segment used to copy
        num_cycles (int): Amount of cycles of PCR to complete
        fall_off_noise (int): Amount of noise to add to generated fall off rate

    Returns:
        List[Tuple[Sequence, Sequence]]: Products of PCR
    """
    products: List[Tuple[Sequence, Sequence]] = [segment]

    for _ in range(num_cycles):
        start_time = time()

        single_strands = denature(products)

        # NOTE May want to calculate fall off rate every elongation
        # instead of every cycle
        fall_off_rate = generate_fall_off_rate(
            d=distance_between_primers(segment, primers), e=fall_off_noise)
        products = annealing_elongation(single_strands, primers, fall_off_rate)
        products = clean_empty_sequences(products)
        products = clean_empty_tuples(products)

        end_time = time()
        time_elapsed = round(end_time - start_time, 2)
        print("Completed cycle {} in {} seconds: {} sequences".format(
            _, time_elapsed, len(products)))

    return products


def denature(sequences: List[Tuple[Sequence, Sequence]]) -> List[Sequence]:
    """Denature a list of paired sequences into single strands

    Args:
        sequences (List[Tuple[Sequence, Sequence]]): List of sequences in PCR reaction

    Returns:
        List[Sequence]: List of single strands
    """
    single_strands = []
    for sequence in sequences:
        if sequence[0] != None:
            single_strands.append(sequence[0])
        if sequence[1] != None:
            single_strands.append(sequence[1])
    return single_strands


def annealing_elongation(single_strands: List[Sequence], primers: Tuple[Sequence, Sequence], fall_off_rate: int) -> List[Sequence]:
    """Anneals and elongates a list of single strands using provided forward and reverse primers

    Args:
        single_strands (List[Sequence]): List of single strands
        primers (Tuple[Sequence, Sequence]): Forward and reverse primers

    Returns:
        List[Sequence]: List of paired sequences
    """

    products: List[Tuple[Sequence, Sequence]] = []

    [f_primer, r_primer] = primers

    for strand in single_strands:
        if strand.bases.find(f_primer.compliment().bases) != -1:
            product = elongate(strand, f_primer, fall_off_rate)
            products.append(product)
            continue

        if strand.bases.find(r_primer.compliment().bases) != -1:
            product = elongate(strand, r_primer, fall_off_rate)
            products.append(product)
            continue

        products.append((strand, None))
    return products


def elongate(template_strand, primer, fall_off_rate) -> Tuple[Sequence, Sequence]:
    """Add nucleotides to single strand to create a paired strand

    Args:
        template_strand ([type]): Single strand to add nucleotides to
        primer ([type]): Location to start adding to
        fall_off_rate ([type]): How many nucleotides should be added starting at the primer

    Returns:
        Tuple[Sequence, Sequence]: Elongated paired strand
    """
    c_primer = primer.compliment()
    index_of_primer = template_strand.bases.find(c_primer.bases)
    primer_length = len(primer.bases)

    start_index = index_of_primer + primer_length

    end_index = max(0, index_of_primer - fall_off_rate)

    copied_segment_sequence = template_strand.bases[end_index:start_index]
    copied_segment = Sequence(copied_segment_sequence).compliment()
    return (template_strand, copied_segment)


def getStats(results: List[Tuple[Sequence, Sequence]]) -> List[int]:
    """
    Find length stats of fragments

    Args: 
        results (List[Tuple[Sequence, Sequence]]): Sequences and their compliments

    Returns:
        fragment_lengths (List[int]): Lengths of each denatured strand
    """

    fragment_lengths: List[int] = []

    M_gene_len = len(results[0][0].bases)

    results.pop(0)

    for fragment in denature(results):
        fragment_lengths.append(len(fragment.bases))

    fragment_lengths.remove(M_gene_len)

    print("Number of Fragments:", len(fragment_lengths))
    print("Max Length:", max(fragment_lengths))
    print("Min Length:", min(fragment_lengths))
    print("Average Length:", (sum(fragment_lengths)/len(fragment_lengths)))

    return fragment_lengths


def graphFragments(fragment_lengths: List[int]):
    """
    Display column graph to visualize distribution and length of fragments

    Args:
        fragment_lengths (List[int]): List containing the lengths of all protein fragments
    """

    i = min(fragment_lengths)
    occurances = []
    length = []
    length_count = []

    # Counts occurances of values and stores in a list
    while i <= max(fragment_lengths):
        occurances.append([i, fragment_lengths.count(i)])
        i += 1

    # Splits list for x, y values
    for item in occurances:
        length.append(item[0])
        length_count.append(item[1])

    plt.bar(length, length_count)
    plt.title("Length and Quantity of M Protein Fragments")
    plt.xlabel("Length of Fragments")
    plt.ylabel("Number of Fragments")
    plt.tight_layout()
    plt.show()
