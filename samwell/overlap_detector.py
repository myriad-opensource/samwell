"""
Utility Classes for Querying Overlaps with Genomic Regions
----------------------------------------------------------

Examples of Detecting Overlaps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    >>> from samwell.overlap_detector import Interval, OverlapDetector
    >>> detector = OverlapDetector()
    >>> query = Interval("chr1", 2, 20)
    >>> detector.overlaps_any(query)
    False
    >>> detector.add(Interval("chr2", 1, 100))
    >>> detector.add(Interval("chr1", 21, 100))
    >>> detector.overlaps_any(query)
    False
    >>> detector.add(Interval("chr1", 1, 1))
    >>> detector.overlaps_any(query)
    True
    >>> detector.get_overlaps(query)
    [Interval("chr1", 1, 1)]
    >>> detector.add(Interval("chr1", 3, 10))
    >>> detector.overlaps_any(query)
    True
    >>> detector.get_overlaps(query)
    [Interval("chr1", 1, 1), interval("chr1", 3, 10)]

Module Contents
~~~~~~~~~~~~~~~

The module contains the following public classes:

    - :class:`~samwell.overlap_detector.Interval` -- Represents a region mapping to the genome
        that is 0-based and open-ended
    - :class:`~samwell.overlap_detector.OverlapDetector` -- Detects and returns overlaps between
        a set of genomic regions and another genomic region
"""

from pathlib import Path
from typing import Dict
from typing import Iterable
from typing import List
from typing import Optional

import attr
from intervaltree import IntervalTree
from pybedtools import BedTool


@attr.s(frozen=True, auto_attribs=True)
class Interval:
    """A region mapping to the genome that is 0-based and open-ended

    Attributes:
        refname (str): the refname (or chromosome)
        start (int): the 0-based start position
        end (int): the 0-based end position (exclusive)
        negative (bool): true if the interval is on the negative strand, false otherwise
        name (Optional[str]): an optional name assigned to the interval
    """
    refname: str = attr.ib()
    start: int = attr.ib()
    end: int = attr.ib()
    negative: bool = False
    name: Optional[str] = None

    def __attrs_post_init__(self) -> None:
        """ Performs simple validation.

        Checks:
            - 0 <= start
            - start < end
        """
        if self.start < 0:
            raise ValueError(f"start is out of range: {self.start}")
        if self.end <= self.start:
            raise ValueError(f"end < start: {self.end} < {self.start}")

    def overlap(self, other: "Interval") -> int:
        """Returns the overlap between this interval and the other, or zero if there is none.

        Args:
            other (Interval): the other interval to find the overlap with
        """
        if self.refname != other.refname:
            return 0

        overlap = min(self.end, other.end) - max(self.start, other.start)
        return overlap if overlap > 0 else 0

    def length(self) -> int:
        """Returns the length of the interval."""
        return self.end - self.start


class OverlapDetector:
    """Detects and returns overlaps between a set of genomic regions and another genomic region.

    If two intervals have the same coordinates, only the first that was added will be kept.

    Since :class:`~samwell.overlap_detector.Interval` objects are used both to populate the
    overlap detector and to query it, the coordinate system in use is also 0-based open-ended.
    """

    def __init__(self) -> None:
        self._refname_to_tree: Dict[str, IntervalTree] = {}

    def add(self, interval: Interval) -> None:
        """Adds a interval to this detector.

        Args:
            interval: the interval to add to this detector
        """
        refname = interval.refname
        if refname not in self._refname_to_tree:
            self._refname_to_tree[refname] = IntervalTree()
        tree = self._get_tree(refname)
        tree[interval.start:interval.end] = interval

    def add_all(self, intervals: Iterable[Interval]) -> None:
        """Adds one or more intervals to this detector.

        Args:
            intervals: the intervals to add to this detector
        """
        for interval in intervals:
            self.add(interval)

    def overlaps_any(self, interval: Interval) -> bool:
        """Determines whether the given interval overlaps any interval in this detector.

           Args:
               interval: the interval to check

           Returns:
               True if and only if the given interval overlaps with any interval in this
               detector.
           """
        tree = self._get_tree(interval.refname)
        if tree is None:
            return False
        else:
            return tree.overlaps(interval.start, interval.end)

    def get_overlaps(self, interval: Interval) -> List[Interval]:
        """Returns any intervals in this detector that overlap the given interval.

          Args:
              interval: the interval to check

          Returns:
              The list of intervals in this detector that overlap the given interval, or the empty
              list if no overlaps exist.  The intervals will be return in ascending genomic order.
          """
        tree = self._get_tree(interval.refname)
        if tree is None:
            return []
        else:
            intervals = [i.data for i in tree.overlap(interval.start, interval.end)]
            return sorted(intervals, key=lambda intv: (intv.start, intv.end))

    def get_enclosed(self, interval: Interval) -> List[Interval]:
        """Returns the set of intervals in this detector that are enclosed by the query
        interval.  I.e. target.start >= query.start and target.end <= query.end.

          Args:
              interval: the query interval

          Returns:
              The list of intervals in this detector that are enclosed within the query interval.
              The intervals will be return in ascending genomic order.
        """
        results = self.get_overlaps(interval)
        return [i for i in results if i.start >= interval.start and i.end <= interval.end]

    def _get_tree(self, refname: str) -> Optional[IntervalTree]:
        """Gets the detector for the given refname

        Args:
            refname: the refname

        Returns:
            the detector for the given refname, or None if there are no regions for that refname
        """
        return self._refname_to_tree.get(refname)

    @classmethod
    def from_bed(cls, path: Path) -> 'OverlapDetector':
        """Builds an :class:`~samwell.overlap_detector.OverlapDetector` from a BED file.

        Args:
            path: the path to the BED file

        Returns:
            An overlap detector for the regions in the BED file.
        """
        detector = OverlapDetector()
        for region in BedTool(str(path)):
            locatable = Interval(refname=region.chrom,
                                 start=region.start,
                                 end=region.end,
                                 negative=region.strand == "-",
                                 name=region.name if region.name != "." else None
                                 )
            detector.add(locatable)
        return detector
