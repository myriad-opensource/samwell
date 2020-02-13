"""Tests for :py:mod:`~samwell.dnautils`"""

import pytest

from samwell import dnautils


def test_reverse_complement() -> None:
    assert dnautils.reverse_complement("") == ""
    assert dnautils.reverse_complement("AATTCCGGaattccgg") == "ccggaattCCGGAATT"
    assert dnautils.reverse_complement("ACGTN") == "NACGT"

    with pytest.raises(KeyError):
        dnautils.reverse_complement("ACGT.GAT")

    with pytest.raises(KeyError):
        dnautils.reverse_complement("RMNY")


def test_mask_long_homopolymers() -> None:
    assert dnautils.mask_long_homopolymers("A", 0) == "N"
    assert dnautils.mask_long_homopolymers("A", 1) == "N"
    assert dnautils.mask_long_homopolymers("A", 2) == "A"
    assert dnautils.mask_long_homopolymers("ACCGGGTTTTAAAAATTTT", 1) == "NNNNNNNNNNNNNNNNNNN"
    assert dnautils.mask_long_homopolymers("ACCGGGTTTTAAAAATTTT", 2) == "ANNNNNNNNNNNNNNNNNN"
    assert dnautils.mask_long_homopolymers("ACCGGGTTTTAAAAATTTT", 3) == "ACCNNNNNNNNNNNNNNNN"
    assert dnautils.mask_long_homopolymers("ACCGGGTTTTAAAAATTTT", 4) == "ACCGGGNNNNNNNNNNNNN"
    assert dnautils.mask_long_homopolymers("ACCGGGTTTTAAAAATTTT", 5) == "ACCGGGTTTTNNNNNTTTT"
    assert dnautils.mask_long_homopolymers("ACCGGGTTTTAAAAATTTT", 6) == "ACCGGGTTTTAAAAATTTT"


def test_has_long_homopolymer() -> None:
    assert dnautils.has_long_homopolymer("A", 0)
    assert not dnautils.has_long_homopolymer("A", 1)
    assert dnautils.has_long_homopolymer("ACCGGGTTTTAAAAATTTT", 4)
    assert not dnautils.has_long_homopolymer("ACCGGGTTTTAAAAATTTT", 5)
    assert not dnautils.has_long_homopolymer("ACCGGGTTTTAAAAATTTT", 10)
