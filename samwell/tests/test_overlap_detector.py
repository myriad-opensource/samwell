"""Tests for :py:mod:`~samwell.overlap_detector`"""

from typing import List

from samwell.overlap_detector import Interval
from samwell.overlap_detector import OverlapDetector


def run_test(targets: List[Interval], query: Interval, results: List[Interval]) -> None:
    detector = OverlapDetector()
    # Use add_all() to covert itself and add()
    detector.add_all(intervals=targets)
    # Test overlaps_any()
    assert detector.overlaps_any(query) == (len(results) > 0)
    # Test get_overlaps()
    assert detector.get_overlaps(query) == results


def test_same_interval() -> None:
    interval = Interval("1", 10, 100)
    run_test(targets=[interval], query=interval, results=[interval])


def test_query_wholly_contained_in_target() -> None:
    target = Interval("1", 10, 100)
    query = Interval("1", 11, 99)
    run_test(targets=[target], query=query, results=[target])


def test_target_wholly_contained_in_query() -> None:
    target = Interval("1", 10, 100)
    query = Interval("1", 9, 101)
    run_test(targets=[target], query=query, results=[target])


def test_target_overlaps_first_base_of_query() -> None:
    target = Interval("1", 10, 100)
    query = Interval("1", 99, 100)
    run_test(targets=[target], query=query, results=[target])


def test_target_overlaps_last_base_of_query() -> None:
    target = Interval("1", 10, 100)
    query = Interval("1", 10, 11)
    run_test(targets=[target], query=query, results=[target])


def test_query_before_target() -> None:
    target = Interval("1", 10, 100)
    query = Interval("1", 9, 10)
    run_test(targets=[target], query=query, results=[])


def test_query_after_target() -> None:
    target = Interval("1", 10, 100)
    query = Interval("1", 100, 101)
    run_test(targets=[target], query=query, results=[])


def test_different_references() -> None:
    target = Interval("1", 10, 100)
    query = Interval("2", 10, 100)
    run_test(targets=[target], query=query, results=[])


def test_multiple_overlaps() -> None:
    interval_a = Interval("1", 10, 20)
    interval_b = Interval("1", 15, 25)
    interval_c = Interval("1", 19, 30)
    interval_d = Interval("1", 24, 35)

    # B overlaps both A and C
    run_test(targets=[interval_a, interval_c], query=interval_b, results=[interval_a, interval_c])
    # C overlaps both A and B
    run_test(targets=[interval_a, interval_b], query=interval_c, results=[interval_a, interval_b])
    # D overlaps only B and C (is after A)
    run_test(targets=[interval_a, interval_b, interval_c],
             query=interval_d,
             results=[interval_b, interval_c])


def test_multiple_references() -> None:
    target_chr1 = Interval("1", 10, 20)
    target_chr2 = Interval("2", 10, 20)
    run_test(targets=[target_chr1, target_chr2], query=target_chr1, results=[target_chr1])
    run_test(targets=[target_chr1, target_chr2], query=target_chr2, results=[target_chr2])


def test_same_interval_twice() -> None:
    interval = Interval("1", 10, 100)
    run_test(targets=[interval, interval], query=interval, results=[interval])


def test_wholly_contained_target() -> None:
    target_inner = Interval("1", 50, 60)
    target_outer = Interval("1", 40, 80)

    run_test(targets=[target_inner, target_outer],
             query=target_inner,
             results=[target_outer, target_inner])


def test_get_enclosed() -> None:
    a = Interval("1", 10, 100)
    b = Interval("1", 15, 20)
    c = Interval("1", 18, 19)
    d = Interval("1", 50, 99)

    detector = OverlapDetector()
    detector.add_all([a, b, c, d])

    assert detector.get_enclosed(Interval("1", 1, 250)) == [a, b, c, d]
    assert detector.get_enclosed(Interval("1", 5, 30)) == [b, c]
    assert detector.get_enclosed(Interval("1", 16, 20)) == [c]
    assert detector.get_enclosed(Interval("1", 15, 19)) == [c]
    assert detector.get_enclosed(Interval("1", 10, 99)) == [b, c, d]
