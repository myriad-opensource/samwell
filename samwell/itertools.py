"""
Functions for Creating Useful Iterators
---------------------------------------

This module contains classes and functions for creating useful iterators.

Examples of a "Peekable" Iterator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"Peekable" iterators are useful to "peek" at the next item in an iterator without consuming it.
For example, this is useful when consuming items in iterator while a predicate is true, and not
consuming the first element where the element is not true.  See the
:func:`~samwell.itertools.PeekableIterator.takewhile` and
:func:`~samwell.itertools.PeekableIterator.dropwhile` methods.

An empty peekable iterator throws StopIteration:

.. code-block:: python

    >>> from samwell.itertools import peekable
    >>> piter = peekable(iter([]))
    >>> piter.peek()
    StopIteration

A peekable iterator will return the next item before consuming it.

.. code-block:: python

    >>> piter = peekable(iter([1, 2, 3]))
    >>> piter.peek()
    1
    >>> next(piter)
    1
    >>> [j for j in piter]
    [2, 3]

The `can_peek()` function can be used to determine if the iterator can be peeked without
StopIteration being thrown:

    >>> piter = peekable([1])
    >>> piter.peek() if piter.can_peek() else -1
    1
    >>> next(piter)
    1
    >>> piter.peek() if piter.can_peek() else -1
    -1
    >>> next(piter)
    StopIteration

The `peekable()` function should be preferred to calling `PeekableIterator`'s constructor
directly as it supports creation from iterable objects as well as iterators, while the constructor
requires an iterator.

Examples of a "Merging" Iterator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A "merging" iterator can merge two iterators in order based on a given ordering function.  This is
useful for merging two iterators that are already in order.

.. code-block:: python

    >>> from samwell.itertools import MergingIterator
    >>> even = iter([2, 4, 6, 8])
    >>> odd = iter([1, 3, 5, 9])
    >>> merging = MergingIterator(even, odd, lambda x: x)
    >>> list(merging)
    [1, 2, 3, 4, 5, 6, 7, 8, 9]

Module Contents
~~~~~~~~~~~~~~~

The module contains the following public classes:

    - :class:`~samwell.itertools.PeekableIterator` -- Iterator that allows you to peek at the
        next value before calling next

    - :class:`samwell.itertools.MergingIterator` -- Iterator that allows merging of two
        iterator using a keyfunc to decide from which iterator to draw the next item

The module contains the following methods:

    - :func:`~samwell.itertools.peekable` -- Creates an iterator that allows you to peek at
        the next value before calling next
"""

from typing import Any
from typing import Callable
from typing import Generic
from typing import Iterable
from typing import Iterator
from typing import List
from typing import TypeVar
from typing import Union


IterType = TypeVar('IterType')


class PeekableIterator(Generic[IterType], Iterator[IterType]):
    """A peekable iterator wrapping an iterable.

    This allows returning the next item without consuming it.

    Args:
        source: an iterator over the objects
    """

    def __init__(self, source: Iterator[IterType]) -> None:
        self._iter: Iterator[IterType] = source
        self._sentinel: Any = object()
        self.__update_peek()

    def __iter__(self) -> Iterator[IterType]:
        return self

    def __next__(self) -> IterType:
        to_return = self.peek()
        self.__update_peek()
        return to_return

    def __update_peek(self) -> None:
        self._peek = next(self._iter, self._sentinel)

    def can_peek(self) -> bool:
        """Returns true if there is a value that can be peeked at, false otherwise."""
        return self._peek is not self._sentinel

    def peek(self) -> IterType:
        """Returns the next element without consuming it, or StopIteration otherwise."""
        if self.can_peek():
            return self._peek
        else:
            raise StopIteration

    def takewhile(self, pred: Callable[[IterType], bool]) -> List[IterType]:
        """Consumes from the iterator while pred is true, and returns the result as a List.

        The iterator is left pointing at the first non-matching item, or if all items match
        then the iterator will be exhausted.

        Args:
            pred: a function that takes the next value from the iterator and returns
                  true or false.

        Returns:
            List[V]: A list of the values from the iterator, in order, up until and excluding
            the first value that does not match the predicate.
        """
        xs: List[IterType] = []
        while self.can_peek() and pred(self._peek):
            xs.append(next(self))
        return xs

    def dropwhile(self, pred: Callable[[IterType], bool]) -> "PeekableIterator[IterType]":
        """Drops elements from the iterator while the predicate is true.

        Updates the iterator to point at the first non-matching element, or exhausts the
        iterator if all elements match the predicate.

        Args:
            pred (Callable[[V], bool]): a function that takes a value from the iterator
            and returns true or false.

        Returns:
            PeekableIterator[V]: a reference to this iterator, so calls can be chained
        """
        while self.can_peek() and pred(self._peek):
            self.__update_peek()
        return self


def peekable(source: Union[Iterator[IterType], Iterable[IterType]]) -> PeekableIterator[IterType]:
    """Creates a peekable iterator that allows you to peek at the next value before calling next

    The peek method will return the next element without consuming it, otherwise StopIteration.

    Args:
        source: either an iterator over the objects, or a callable that is called until it
            returns the sentinel.

    Returns:
        a :class:`~samwell.itertools.PeekableIterator`
    """
    return PeekableIterator(source=iter(source))


class MergingIterator(Generic[IterType], Iterator[IterType]):
    """An iterator that merges two iterators; if they are sorted and keyfunc is passed, yields
    results in order.

    Args:
        iter1: an iterator
        iter2: an iterator
        keyfunc: a function that extracts a key from an item that is used to order items
    """

    def __init__(self,
                 iter1: Iterator[IterType],
                 iter2: Iterator[IterType],
                 keyfunc: Callable[[IterType], Any]) -> None:
        self._iter1 = peekable(iter1)
        self._iter2 = peekable(iter2)
        self._keyfunc = keyfunc

    def __iter__(self) -> Iterator[IterType]:
        return self

    def __next__(self) -> IterType:
        if self._iter1.can_peek() and self._iter2.can_peek():
            k1 = self._keyfunc(self._iter1.peek())
            k2 = self._keyfunc(self._iter2.peek())
            return next(self._iter1 if k1 <= k2 else self._iter2)
        elif self._iter1.can_peek():
            return next(self._iter1)
        elif self._iter2.can_peek():
            return next(self._iter2)
        else:
            raise StopIteration
