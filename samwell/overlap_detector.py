"""
Utility Classes for Querying Overlaps with Genomic Regions
----------------------------------------------------------

DEPRECATED - if you have the option use `~pybedlite.overlap_detector` in favor of this.

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


from pybedlite.overlap_detector import Interval
from pybedlite.overlap_detector import OverlapDetector

__all__ = ["Interval", "OverlapDetector"]
