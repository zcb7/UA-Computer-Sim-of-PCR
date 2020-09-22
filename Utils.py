from random import randrange
from typing import Tuple

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
