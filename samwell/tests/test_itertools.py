"""Tests for :py:mod:`~samwell.itertools`"""

import pytest

from samwell.itertools import PeekableIterator
from samwell.itertools import peekable
from samwell.itertools import MergingIterator


def test_peekable_iterator_empty() -> None:
    empty_iter: PeekableIterator[None] = peekable([])
    assert not empty_iter.can_peek()
    assert empty_iter.maybe_peek() is None, "maybe_peek was not None for empty iterator"
    with pytest.raises(StopIteration):
        empty_iter.peek()
    with pytest.raises(StopIteration):
        next(empty_iter)


def test_peekable_iterator_nonempty() -> None:
    nonempty_iter = peekable(range(10))
    for i in range(10):
        assert nonempty_iter.can_peek()
        assert nonempty_iter.peek() == i
        assert nonempty_iter.maybe_peek() == i, "maybe_peek value didn't match expectation"
        assert next(nonempty_iter) == i

    assert nonempty_iter.maybe_peek() is None, "maybe_peek was not None for exhausted iterator"
    with pytest.raises(StopIteration):
        nonempty_iter.peek()
    with pytest.raises(StopIteration):
        next(nonempty_iter)


def test_peekable_with_nones() -> None:
    xs = [1, 2, None, 4, None, 6]
    iterator = peekable(xs)

    for i in range(len(xs)):
        assert iterator.peek() is xs[i]
        assert iterator.maybe_peek() is xs[i]
        assert next(iterator) is xs[i]


def test_takewhile() -> None:
    xs = [2, 4, 6, 8, 11, 13, 15, 17, 19, 20, 22, 24]
    iterator = peekable(xs)
    assert iterator.takewhile(lambda x: x % 2 == 0) == [2, 4, 6, 8]
    assert iterator.takewhile(lambda x: x % 2 == 1) == [11, 13, 15, 17, 19]
    assert iterator.takewhile(lambda x: x % 2 == 1) == []
    assert iterator.takewhile(lambda x: x % 2 == 0) == [20, 22, 24]


def test_dropwhile() -> None:
    xs = [2, 4, 6, 8, 11, 13, 15, 17, 19, 20, 22, 24]
    iterator = peekable(xs)
    iterator.dropwhile(lambda x: x % 2 == 0)
    iterator.dropwhile(lambda x: x <= 20)
    assert list(iterator) == [22, 24]


def test_merging_iterator() -> None:
    xs = [1, 3, 5, 7, 9]
    ys = [2, 4, 6, 8, 9]
    ms = MergingIterator(iter(xs), iter(ys), keyfunc=lambda x: x)
    assert list(ms) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 9]

    s1 = ["one", "enormous", "hippopotamus"]
    s2 = ["a", "little", "diplodocus"]
    ss = MergingIterator(iter(s1), iter(s2), keyfunc=lambda x: len(x))
    assert list(ss) == ["a", "one", "little", "enormous", "diplodocus", "hippopotamus"]
