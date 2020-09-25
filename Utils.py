from random import randrange
from typing import List, Tuple

from GenBank import Sequence


def generate_fall_off_rate(d=200, e=50):
    """Generate a random integer fall-off rate for Taq Polymerase

    Args:
        d (int, optional): The constant pivot. Defaults to 200.
        e (int, optional): The maximum / minimum random noise. Must be positive and non-zero. Defaults to 50.

    Returns:
        [type]: A random fall off rate
    """
    return randrange(d - e, d + e)


def distance_between_primers(segment: Tuple[Sequence, Sequence], primers: Tuple[Sequence, Sequence]):
    """Computes the distance between forward and reverse primers on a segment

    Args:
        segment (Tuple[Sequence, Sequence]): The segment that both primers bind to
        primers (Tuple[Sequence, Sequence]): The forward and reverse primers

    Returns:
        [type]: Distance between forward and reverse primers
    """
    start_index = segment[0].bases.find(primers[0].bases)
    end_index = len(segment[1].bases) - segment[1].bases.find(primers[1].bases)
    return end_index - start_index


def clean_empty_sequences(sequence_pairs: List[Tuple[Sequence, Sequence]]) -> List[Tuple]:
    """Replace empty sequences with None type

    Args:
        sequence_pairs (List[Tuple[Sequence, Sequence]]): A list of sequence pairs

    Returns:
        List[Tuple]: A cleaned list of sequence pairs
    """
    for index, pair in enumerate(sequence_pairs):
        template = pair[0]
        compliment = pair[1]

        if pair[0] is not None and len(pair[0].bases) == 0:
            template = None
        if pair[1] is not None and len(pair[1].bases) == 0:
            compliment = None

        sequence_pairs[index] = (template, compliment)

    return sequence_pairs


def clean_empty_tuples(tuples: List[Tuple]) -> List[Tuple]:
    """Remove tuples of consisting of (None, None)

    Args:
        tuples (List[Tuple]): List of sequence pairs

    Returns:
        List[Tuple]: A cleaned list of sequence pairs
    """
    filtered_tuples = filter(lambda x: not (
        x[0] == None and x[1] == None), tuples)
    return list(filtered_tuples)
