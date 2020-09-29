from time import time
from typing import List, Tuple

from matplotlib import pyplot as plt

from GenBank import Sequence
from Utils import (clean_empty_sequences, clean_empty_tuples,
                   distance_between_primers, generate_fall_off_rate)


def PCR(segment: Tuple[Sequence, Sequence], primers: Tuple[Sequence, Sequence], num_cycles: int = 20, fall_off_rate_pivot: int = None, fall_off_noise: int = None) -> List[Tuple[Sequence, Sequence]]:
    """Simulate Polymerase Chain Reaction

    Args:
        segment (Tuple[Sequence, Sequence]): Segment to copy
        primers (Tuple[Sequence, Sequence]): Primers on segment used to copy
        num_cycles (int): Amount of cycles of PCR to complete
        fall_off_rate_pivot (int): The pivot used to generate the fall off rate
        fall_off_noise (int): Amount of noise to add to generated fall off rate

    Returns:
        List[Tuple[Sequence, Sequence]]: Products of PCR
    """
    products: List[Tuple[Sequence, Sequence]] = [segment]

    for _ in range(num_cycles):
        start_time = time()

        single_strands = denature(products)

        products = annealing_elongation(
            single_strands, primers, fall_off_rate_pivot, fall_off_noise)
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


def annealing_elongation(single_strands: List[Sequence], primers: Tuple[Sequence, Sequence], fall_off_rate_pivot: int, fall_off_noise: int) -> List[Sequence]:
    """Anneals and elongates a list of single strands using provided forward and reverse primers

    Args:
        single_strands (List[Sequence]): List of single strands
        primers (Tuple[Sequence, Sequence]): Forward and reverse primers
        fall_off_rate_pivot (int): The pivot used to generate the fall off rate
        fall_off_noise (int): The amount of noise to apply to the fall off rate

    Returns:
        List[Sequence]: List of paired sequences
    """

    products: List[Tuple[Sequence, Sequence]] = []

    [f_primer, r_primer] = primers

    for strand in single_strands:
        if strand.bases.find(f_primer.compliment().bases) != -1:
            product = elongate(
                strand, f_primer, fall_off_rate_pivot, fall_off_noise)
            products.append(product)
            continue

        if strand.bases.find(r_primer.compliment().bases) != -1:
            product = elongate(
                strand, r_primer, fall_off_rate_pivot, fall_off_noise)
            products.append(product)
            continue

        products.append((strand, None))
    return products


def elongate(template_strand, primer, fall_off_rate_pivot, fall_off_noise) -> Tuple[Sequence, Sequence]:
    """Add nucleotides to single strand to create a paired strand

    Args:
        template_strand ([type]): Single strand to add nucleotides to
        primer ([type]): Location to start adding to
        fall_off_rate_pivot (int): The pivot used to generate the fall off rate
        fall_off_rate ([type]): The amount of noise to apply to the fall off rate

    Returns:
        Tuple[Sequence, Sequence]: Elongated paired strand
    """
    c_primer = primer.compliment()
    index_of_primer = template_strand.bases.find(c_primer.bases)
    primer_length = len(primer.bases)

    start_index = index_of_primer + primer_length

    fall_off_rate = generate_fall_off_rate(
        d=fall_off_rate_pivot, e=fall_off_noise)
    end_index = max(0, index_of_primer - fall_off_rate)

    copied_segment_sequence = template_strand.bases[end_index:start_index]
    copied_segment = Sequence(copied_segment_sequence).compliment()
    return (template_strand, copied_segment)


def getStats(results: List[Tuple[Sequence, Sequence]], primers: Tuple[Sequence, Sequence]) -> List[int]:
    """
    Find length stats of fragments

    Args: 
        results (List[Tuple[Sequence, Sequence]]): Sequences and their compliments

    Returns:
        fragment_lengths (List[int]): Lengths of each denatured strand
    """

    fragment_lengths: List[int] = []

    original_strand = results[0][0]
    M_gene_len = len(original_strand.bases)

    for fragment in denature(results):
        fragment_lengths.append(len(fragment.bases))

    fragment_lengths = list(
        filter(lambda x: x != M_gene_len, fragment_lengths))

    copied_segment_length = distance_between_primers(
        (original_strand, original_strand.compliment()), primers)
    copied_segment_start = original_strand.bases.index(primers[0].bases)
    copied_segment = Sequence(original_strand.bases[copied_segment_start:
                                           copied_segment_start + copied_segment_length])

    print("Number of Fragments:", len(fragment_lengths))
    print("Max Length:", max(fragment_lengths))
    print("Min Length:", min(fragment_lengths))
    print("Average Length:", (sum(fragment_lengths)/len(fragment_lengths)))
    print("GC Content for M gene:", original_strand.gc_content())
    print("GC Content for copied segment:", copied_segment.gc_content())
    print("Copied segment length:", copied_segment_length)

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

    plt.bar(length, length_count)
    plt.title("Length and Quantity of M Protein Fragments (Logarithmic)")
    plt.xlabel("Length of Fragments")
    plt.ylabel("Number of Fragments")
    plt.yscale("log")
    plt.tight_layout()
    plt.show()
