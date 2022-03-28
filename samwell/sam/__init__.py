"""
Utility Classes and Methods for SAM/BAM
---------------------------------------

This module contains utility classes for reading and writing SAM/BAM files, as well as for
manipulating Cigars.  It is recommended to use the :func:`~samwell.sam.reader` and
:func:`~samwell.sam.writer` methods rather than :class:`pysam.AlignmentFile` directly (see
below for motivation).

Motivation for Reader and Writer methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following are the reasons for choosing to implement methods to open a SAM/BAM file for
reading and writing, rather than relying on :class:`pysam.AlignmentFile` directly:

1. Provides a centralized place for the implementation of opening a SAM/BAM for reading and
   writing.  This is useful if any additional parameters are added, or changes to standards or
   defaults are made.
2. Makes the requirement to provide a header when opening a file for writing more explicit.
3. Adds support for :class:`~pathlib.Path`.
4. Remove the reliance on specifying the mode correctly, including specifying the file type (i.e.
   SAM, BAM, or CRAM), as well as additional options (ex. compression level).  This makes the
   code more explicit and easier to read.
5. An explicit check is performed to ensure the file type is specified when writing using a
   file-like object rather than a path to a file.

Examples of Opening a SAM/BAM for Reading or Writing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Opening a SAM/BAM file for reading, auto-recognizing the file-type by the file extension.  See
:class:`~samwell.sam.SamFileType` for the supported file types.

.. code-block:: python

    >>> from samwell.sam import reader
    >>> with reader("/path/to/sample.sam") as fh:
    ...     for record in fh:
    ...         print(record.name)  # do something
    >>> with reader("/path/to/sample.bam") as fh:
    ...     for record in fh:
    ...         print(record.name)  # do something

Opening a SAM/BAM file for reading, explicitly passing the file type.

    >>> from samwell.sam import SamFileType
    >>> with reader(path="/path/to/sample.ext1", file_type=SamFileType.SAM) as fh:
    ...     for record in fh:
    ...         print(record.name)  # do something
    >>> with reader(path="/path/to/sample.ext2", file_type=SamFileType.BAM) as fh:
    ...     for record in fh:
    ...         print(record.name)  # do something

Opening a SAM/BAM file for reading, using an existing file-like object

    >>> with open("/path/to/sample.sam", "rb") as file_object:
    ...     with reader(path=file_object, file_type=SamFileType.BAM) as fh:
    ...         for record in fh:
    ...             print(record.name)  # do something

Opening a SAM/BAM file for writing follows similar to the :func:`~samwell.sam.reader` method,
but the SAM file header object is required.

    >>> from samwell.sam import writer
    >>> header: Dict[str, Any] = {
    ...     "HD": {"VN": "1.5", "SO": "coordinate"},
    ...     "RG": [{"ID": "1", "SM": "1_AAAAAA", "LB": "lib", "PL": "ILLUMINA", "PU": "xxx.1"}],
    ...     "SQ":  [
    ...         {"SN": "chr1", "LN": 249250621},
    ...         {"SN": "chr2", "LN": 243199373}
    ...     ]
    ... }
    >>> with writer(path="/path/to/sample.bam", header=header) as fh:
    ...     pass  # do something

Examples of Manipulating Cigars
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Creating a :class:`~samwell.sam.Cigar` from a :class:`pysam.AlignedSegment`.

    >>> from samwell.sam import Cigar
    >>> with reader("/path/to/sample.sam") as fh:
    ...     record = next(fh)
    ...     cigar = Cigar.from_cigartuples(record.cigartuples)
    ...     print(str(cigar))
    50M2D5M10S

Creating a :class:`~samwell.sam.Cigar` from a :class:`str`.

    >>> cigar = Cigar.from_cigarstring("50M2D5M10S")
    >>> print(str(cigar))
    50M2D5M10S

If the cigar string is invalid, the exception message will show you the problem character(s) in
square brackets.

    >>> cigar = Cigar.from_cigarstring("10M5U")
    ... CigarException("Malformed cigar: 10M5[U]")

The cigar contains a tuple of :class:`~samwell.sam.CigarElement`s.  Each element contains the
cigar operator (:class:`~samwell.sam.CigarOp`) and associated operator length.  A number of
useful methods are part of both classes.

The number of bases aligned on the query (i.e. the number of bases consumed by the cigar from
the query):

    >>> cigar = Cigar.from_cigarstring("50M2D5M2I10S")
    >>> [e.length_on_query for e in cigar.elements]
    [50, 0, 5, 2, 10]
    >>> [e.length_on_target for e in cigar.elements]
    [50, 2, 5, 0, 0]
    >>> [e.operator.is_indel for e in cigar.elements]
    [False, True, False, True, False]

Module Contents
~~~~~~~~~~~~~~~

The module contains the following public classes:

    - :class:`~samwell.sam.SamFileType` -- Enumeration of valid SAM/BAM/CRAM file types.
    - :class:`~samwell.sam.SamOrder` -- Enumeration of possible SAM/BAM/CRAM sort orders.
    - :class:`~samwell.sam.CigarOp` -- Enumeration of operators that can appear in a Cigar string.
    - :class:`~samwell.sam.CigarElement` -- Class representing an element in a Cigar string.
    - :class:`~samwell.sam.CigarParsingException` -- The exception raised specific to parsing a
        cigar
    - :class:`~samwell.sam.Cigar` -- Class representing a cigar string.

The module contains the following methods:

    - :func:`~samwell.sam.reader` -- opens a SAM/BAM/CRAM file for reading.
    - :func:`~samwell.sam.writer` -- opens a SAM/BAM/CRAM file for writing
    - :func:`~samwell.sam.set_qc_fail` -- sets the QC fail flag in a
        :class:`pysam.AlignedSegment` record and sets additional SAM tags giving the tool name and
        reason for why the QC fail flag was set.
    - :func:`~samwell.sam.get_qc_fail` -- gets the tool name and reason for why the QC fail flag
        was set, or None if it is not set.
"""

import enum
import io
from pathlib import Path
from typing import Any
from typing import Callable
from typing import Dict
from typing import IO
from typing import List
from typing import Optional
from typing import Tuple
from typing import Union
from typing import TYPE_CHECKING
import attr
import pysam
import sys
from pysam import AlignmentFile as SamFile
from pysam import AlignmentHeader as SamHeader
from pysam import AlignedSegment

if TYPE_CHECKING or sys.version_info < (3, 8, 0):
    from typing_extensions import Final
else:
    from typing import Final

"""The valid base classes for opening a SAM/BAM/CRAM file."""
SamPath = Union[IO[Any], Path, str]

"""The reference index to use to indicate no reference in SAM/BAM."""
NO_REF_INDEX: int = -1

"""The reference name to use to indicate no reference in SAM/BAM."""
NO_REF_NAME: str = "*"

"""The reference position to use to indicate no position in SAM/BAM."""
NO_REF_POS: int = -1


@enum.unique
class SamFileType(enum.Enum):
    """Enumeration of valid SAM/BAM/CRAM file types.

    Attributes:
        mode (str): The additional mode character to add when opening this file type.
        ext (str): The standard file extension for this file type.
    """

    def __init__(self, mode: str, ext: str) -> None:
        self.mode: Final[str] = mode
        self.ext: Final[str] = ext

    SAM = ("", ".sam")
    BAM = ("b", ".bam")
    CRAM = ("c", ".cram")

    @classmethod
    def from_path(cls, path: Union[Path, str]) -> 'SamFileType':
        """Infers the file type based on the file extension.

        Args:
            path: the path to the SAM/BAM/CRAM to read or write.
        """
        ext = Path(path).suffix
        try:
            return next(iter([tpe for tpe in SamFileType if tpe.ext == ext]))
        except StopIteration:
            raise ValueError(f"Could not infer file type from {path}")


"""The classes that should be treated as file-like classes"""
_IOClasses = (
    io.TextIOBase,
    io.BufferedIOBase,
    io.RawIOBase,
    io.IOBase
)


def _pysam_open(path: SamPath,
                open_for_reading: bool,
                file_type: Optional[SamFileType] = None,
                **kwargs: Any) -> SamFile:
    """Opens a SAM/BAM/CRAM for reading or writing.

    Args:
        path: a file handle or path to the SAM/BAM/CRAM to read or write.
        open_for_reading: True to open for reading, false otherwise.
        file_type: the file type to assume when opening the file.  If None, then the file type
            will be auto-detected for reading and must be a path-like object for writing.
        kwargs: any keyword arguments to be passed to
        :class:`~pysam.AlignmentFile`; may not include "mode".
    """

    if isinstance(path, (str, Path)):  # type: ignore
        file_type = file_type or SamFileType.from_path(path)
        path = str(path)
    elif not isinstance(path, _IOClasses):  # type: ignore
        open_type = "reading" if open_for_reading else "writing"
        raise TypeError(f"Cannot open '{type(path)}' for {open_type}.")

    if file_type is None and not open_for_reading:
        raise ValueError("file_type must be given when writing to a file-like object")

    # file_type must be set when writing, so if file_type is None, then we must be opening it
    # for reading.  Hence, only set mode in kwargs to pysam when file_type is set and when
    # writing since we can let pysam auto-recognize the file type when reading.  See discussion:
    # https://github.com/pysam-developers/pysam/issues/655
    if file_type is not None:
        kwargs["mode"] = "r" if open_for_reading else "w" + file_type.mode
    else:
        assert open_for_reading, "Bug: file_type was None but open_for_reading was False"

    # Open it!
    return pysam.AlignmentFile(path, **kwargs)


def reader(path: SamPath,
           file_type: Optional[SamFileType] = None
           ) -> SamFile:
    """Opens a SAM/BAM/CRAM for reading.

        Args:
            path: a file handle or path to the SAM/BAM/CRAM to read or write.
            file_type: the file type to assume when opening the file.  If None, then the file
                type will be auto-detected.
       """
    return _pysam_open(path=path, open_for_reading=True, file_type=file_type)


def writer(path: SamPath,
           header: Union[str, Dict[str, Any], SamHeader],
           file_type: Optional[SamFileType] = None) -> SamFile:
    """Opens a SAM/BAM/CRAM for writing.

        Args:
            path: a file handle or path to the SAM/BAM/CRAM to read or write.
            header: Either a string to use for the header or a multi-level dictionary.  The
                multi-level dictionary should be given as follows.  The first level are the four
                types (‘HD’, ‘SQ’, ...). The second level are a list of lines, with each line being
                a list of tag-value pairs. The header is constructed first from all the defined
                fields, followed by user tags in alphabetical order.
            file_type: the file type to assume when opening the file.  If None, then the
                filetype will be auto-detected and must be a path-like object.
          """
    # Set the header for pysam's AlignmentFile
    key = "text" if isinstance(header, str) else "header"
    kwargs = {key: header}

    return _pysam_open(path=path, open_for_reading=False, file_type=file_type, **kwargs)


class _CigarOpUtil:
    """Some useful constants to speed up methods on CigarOp"""

    """A dictionary from the cigar op code to the cigar op char.

    This is to speed up the translation of cigar op code to CigarOp in CigarOp, so needs to be
    declared beforehand.
    """
    CODE_TO_CHARACTER: Dict[int, str] = {0: "M", 1: "I", 2: "D", 3: "N", 4: "S", 5: "H", 6: "P",
                                         7: "EQ", 8: "X"}


@enum.unique
class CigarOp(enum.Enum):
    """Enumeration of operators that can appear in a Cigar string.

    Attributes:
        code (int): The :py:mod:`~pysam` cigar operator code.
        character (int): The single character cigar operator.
        consumes_query (bool): True if this operator consumes query bases, False otherwise.
        consumes_target (bool): True if this operator consumes target bases, False otherwise.
    """

    M = (0, 'M', True, True)  #: Match or Mismatch the reference
    I = (1, 'I', True, False)  #: Insertion versus the reference  # noqa: E741
    D = (2, 'D', False, True)  #: Deletion versus the reference
    N = (3, 'N', False, True)  #: Skipped region from the reference
    S = (4, 'S', True, False)  #: Soft clip
    H = (5, 'H', False, False)  #: Hard clip
    P = (6, 'P', False, False)  #: Padding
    EQ = (7, '=', True, True)  #: Matches the reference
    X = (8, 'X', True, True)  #: Mismatches the reference

    def __init__(self,
                 code: int,
                 character: str,
                 consumes_query: bool,
                 consumes_reference: bool) -> None:
        self.code = code
        self.character = character
        self.consumes_query = consumes_query
        self.consumes_reference = consumes_reference

    @staticmethod
    def from_character(character: str) -> 'CigarOp':
        """Returns the operator from the single character."""
        if CigarOp.EQ.character == character:
            return CigarOp.EQ
        else:
            return CigarOp[character]

    @staticmethod
    def from_code(code: int) -> 'CigarOp':
        """Returns the operator from the given operator code.

        Note: this is mainly used to get the operator from :py:mod:`~pysam`.
        """
        return CigarOp[_CigarOpUtil.CODE_TO_CHARACTER[code]]

    @property
    def is_indel(self) -> bool:
        """Returns true if the operator is an indel, false otherwise. """
        return self == CigarOp.I or self == CigarOp.D


@attr.s(frozen=True, slots=True)
class CigarElement:
    """ Represents an element in a Cigar

    Attributes:
        - length (int): the length of the element
        - operator (CigarOp): the operator of the element
    """

    length: int = attr.ib()
    operator: CigarOp = attr.ib()

    @length.validator
    def _validate_length(self, attribute: Any, value: int) -> None:
        """Validates the length attribute is greater than zero."""
        if value <= 0:
            raise ValueError(f"Cigar element must have a length > 0, found {value}")

    @property
    def length_on_query(self) -> int:
        """Returns the length of the element on the query sequence."""
        return self.length if self.operator.consumes_query else 0

    @property
    def length_on_target(self) -> int:
        """Returns the length of the element on the target (often reference) sequence."""
        return self.length if self.operator.consumes_reference else 0

    def __str__(self) -> str:
        return f"{self.length}{self.operator.character}"


class CigarParsingException(Exception):
    """The exception raised specific to parsing a cigar."""
    pass


@attr.s(frozen=True, slots=True)
class Cigar:
    """Class representing a cigar string.

    Attributes:
        - elements (Tuple[CigarElement, ...]): zero or more cigar elements
    """

    elements: Tuple[CigarElement, ...] = attr.ib(default=())

    @classmethod
    def from_cigartuples(cls, cigartuples: Optional[List[Tuple[int, int]]]) -> 'Cigar':
        """Returns a Cigar from a list of tuples returned by pysam.

        Each tuple denotes the operation and length.  See
        :class:`~samwell.sam.CigarOp` for more information on the
        various operators.  If None is given, returns an empty Cigar.
        """
        if cigartuples is None or cigartuples == []:
            return Cigar()
        try:
            elements = []
            for code, length in cigartuples:
                operator = CigarOp.from_code(code)
                elements.append(CigarElement(length, operator))
            return Cigar(tuple(elements))
        except Exception as ex:
            raise CigarParsingException(f"Malformed cigar tuples: {cigartuples}") from ex

    @classmethod
    def _pretty_cigarstring_exception(cls,
                                      cigarstring: str,
                                      index: int) -> CigarParsingException:
        """Raises an exception highlighting the malformed character"""
        prefix = cigarstring[:index]
        character = cigarstring[index] if index < len(cigarstring) else ""
        suffix = cigarstring[index + 1:]
        pretty_cigarstring = f"{prefix}[{character}]{suffix}"
        message = f"Malformed cigar: {pretty_cigarstring}"
        return CigarParsingException(message)

    @classmethod
    def from_cigarstring(cls, cigarstring: str) -> 'Cigar':
        """Constructs a Cigar from a string returned by pysam.

        If "*" is given, returns an empty Cigar.
        """
        if cigarstring == "*":
            return Cigar()

        cigarstring_length = len(cigarstring)
        if cigarstring_length == 0:
            raise CigarParsingException("Cigar string was empty")

        elements = []
        i = 0
        while i < cigarstring_length:
            if not cigarstring[i].isdigit():
                raise cls._pretty_cigarstring_exception(cigarstring, i)  # type: ignore
            length = int(cigarstring[i])
            i += 1
            while i < cigarstring_length and cigarstring[i].isdigit():
                length = (length * 10) + int(cigarstring[i])
                i += 1
            if i == cigarstring_length:
                raise cls._pretty_cigarstring_exception(cigarstring, i)  # type: ignore
            try:
                operator = CigarOp.from_character(cigarstring[i])
                elements.append(CigarElement(length, operator))
            except KeyError as ex:
                # cigar operator was not valid
                raise cls._pretty_cigarstring_exception(cigarstring, i) from ex  # type: ignore
            except IndexError as ex:
                # missing cigar operator (i == len(cigarstring))
                raise cls._pretty_cigarstring_exception(cigarstring, i) from ex  # type: ignore
            i += 1
        return Cigar(tuple(elements))

    def __str__(self) -> str:
        if self.elements:
            return "".join([str(e) for e in self.elements])
        else:
            return "*"

    def reversed(self) -> "Cigar":
        """Returns a copy of the Cigar with the elements in reverse order."""
        return Cigar(tuple(reversed(self.elements)))

    def length_on_query(self) -> int:
        """Returns the length of the alignment on the query sequence."""
        return sum([elem.length_on_query for elem in self.elements])

    def length_on_target(self) -> int:
        """Returns the length of the alignment on the target sequence."""
        return sum([elem.length_on_target for elem in self.elements])

    def coalesce(self) -> "Cigar":
        """Returns a copy of the cigar adjacent operators of the same type coalesced into single
        operators."""
        new_elements: List[CigarElement] = []
        element_index: int = 0
        while element_index < len(self.elements):
            cur_element: CigarElement = self.elements[element_index]
            op_length: int = cur_element.length
            element_index += 1
            while (element_index < len(self.elements) and
                    cur_element.operator == self.elements[element_index].operator):
                op_length += self.elements[element_index].length
                element_index += 1
            new_elements.append(CigarElement(operator=cur_element.operator, length=op_length))
        return Cigar(tuple(new_elements))


# The SAM tag to store which tool caused the QC fail flag to be set
QcFailToolTag = 'qt'


# The SAM tag to store the reason why the tool caused the QC flag to be set
QcFailReasonTag = 'qr'


def set_qc_fail(rec: pysam.AlignedSegment, tool: Callable[..., Any], reason: str) -> None:
    """Sets the QC fail flag, and adds tags containing the tool name and reason for failing.
    Args:
        rec: the record to fail
        tool: the tool (as a callable) that failed this record
        reason: the reason for failing
    """
    assert '\t' not in reason, f"Reason may not contain tabs: {reason}"
    rec.is_qcfail = True
    rec.set_tag(QcFailToolTag, tool.__name__)
    rec.set_tag(QcFailReasonTag, reason)


def get_qc_fail(rec: pysam.AlignedSegment) -> Optional[Tuple[str, str]]:
    """Gets the tool and reason for why the QC fail flag is set, otherwise None if not set.

    If the QC fail flag is set, but the tool and filter reason SAM tags are not set, None will be
    returned.  Use pysam.AlignedSegment.is_qcfail() to check if the record is simply QC failed.

    Args:
        rec: the record to fail
    """
    if not rec.is_qcfail or not rec.has_tag(QcFailToolTag):
        return None
    else:
        tool_value = rec.get_tag(QcFailToolTag)
        reason_value = rec.get_tag(QcFailReasonTag)
        return (tool_value, reason_value)


def get_qc_fail_by_tool(rec: pysam.AlignedSegment,
                        tool: Callable[..., Any] = None) -> Optional[Tuple[str, str]]:
    """Gets the tool and reason for why the QC fail flag if the flag was set by the passed tool.

    None will be returned in the following cases:
      - The QC fail flag is not set
      - The QC fail flag isset, but the tool and filter reason SAM tags are not set
      - The tool and filter reason SAM tags were set by a different tool

    Use pysam.AlignedSegment.is_qcfail() to check if the record is simply QC failed.

    Args:
        rec: the record to fail
        tool: the tool that must have set the QC fail flag
    """
    maybe_tool_and_reason = get_qc_fail(rec)
    if maybe_tool_and_reason is None:
        return maybe_tool_and_reason
    else:
        tool_value = maybe_tool_and_reason[0]
        return maybe_tool_and_reason if tool.__name__ == tool_value else None


def isize(r1: AlignedSegment, r2: AlignedSegment) -> int:
    """Computes the insert size for a pair of records."""
    if r1.is_unmapped or r2.is_unmapped or r1.reference_id != r2.reference_id:
        return 0
    else:
        r1_pos = r1.reference_end if r1.is_reverse else r1.reference_start
        r2_pos = r2.reference_end if r2.is_reverse else r2.reference_start
        return r2_pos - r1_pos


def set_pair_info(r1: AlignedSegment, r2: AlignedSegment, proper_pair: bool = True) -> None:
    """Resets mate pair information between reads in a pair. Requires that both r1
    and r2 are mapped.  Can be handed reads that already have pairing flags setup or
    independent R1 and R2 records that are currently flagged as SE reads.

    Args:
        r1: read 1
        r2: read 2 with the same queryname as r1
    """
    assert not r1.is_unmapped, f"Cannot process unmapped mate {r1.query_name}/1"
    assert not r2.is_unmapped, f"Cannot process unmapped mate {r2.query_name}/2"
    assert r1.query_name == r2.query_name, (
        f"Attempting to pair reads with different qnames {r1.query_name} vs {r2.query_name}."
    )

    for r in [r1, r2]:
        r.is_paired = True
        r.is_proper_pair = proper_pair

    r1.is_read1 = True
    r1.is_read2 = False
    r2.is_read2 = True
    r2.is_read1 = False

    for src, dest in [(r1, r2), (r2, r1)]:
        dest.next_reference_id = src.reference_id
        dest.next_reference_start = src.reference_start
        dest.mate_is_reverse = src.is_reverse
        dest.mate_is_unmapped = False
        dest.set_tag("MC", src.cigarstring)

    insert_size = isize(r1, r2)
    r1.template_length = insert_size
    r2.template_length = - insert_size


@enum.unique
class SamOrder(enum.Enum):
    """
    Enumerations of possible sort orders for a SAM file.
    """

    Unsorted = "unsorted"  #: the SAM / BAM / CRAM is unsorted
    Coordinate = "coordinate"  #: coordinate sorted
    QueryName = "queryname"  #: queryname sorted
    Unknown = "unknown"  # Unknown SAM / BAM / CRAM sort order
