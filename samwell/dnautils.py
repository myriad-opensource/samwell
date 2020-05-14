"""
Utility Functions for Manipulating DNA sequences.
-------------------------------------------------

This module contains utility functions for manipulating DNA sequences.

"""

from typing import Dict

_RC_DICT: Dict[str, str] = (
    dict(A='T', C='G', G='C', T='A', a='t', c='g', g='c', t='a', N='N')
)


def reverse_complement(bases: str) -> str:
    """Reverse complements a base sequence.

    Arguments:
        bases: the bases to be reverse complemented.

    Returns:
        the reverse complement of the provided base string
    """
    return ''.join([_RC_DICT[b] for b in bases[::-1]])


def mask_long_homopolymers(bases: str, min_long_hp_length: int, mask_base: str = 'N') -> str:
    """Returns the bases masked for regions with long homopolymers

    Args:
        bases: the bases to mask.
        min_long_hp_length: the minimum homopolymer length (inclusive) to mask.
        mask_base: the base to use when masking
    """
    masked = list(bases)
    count = 1
    last_base = bases[0]
    for i in range(1, len(bases)):
        cur_base = bases[i]
        if last_base == cur_base:
            count += 1
        else:
            if count >= min_long_hp_length:
                masked[i - count:i] = mask_base * count
            last_base = cur_base
            count = 1
    if count >= min_long_hp_length:
        masked[-count:] = mask_base * count
    return ''.join(masked)


def has_long_homopolymer(bases: str, max_hp_length: int) -> bool:
    '''Returns true if the given bases has a homopolymer length longer than the given length.

    Args:
        bases: the bases to examine.
        max_hp_length: the maximum homopolymer length to allow.
    '''
    count = 1
    last_base = bases[0]
    for i in range(1, len(bases)):
        cur_base = bases[i]
        if last_base == cur_base:
            count += 1
            if count > max_hp_length:
                return True
        else:
            last_base = cur_base
            count = 1
    return count > max_hp_length
