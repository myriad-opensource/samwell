"""
Utility methods for running BWA
-------------------------------

This module contains methods for running BWA.  Currently only the "mem" algorithm is supported:
    - :func:`~samwell.bwa_mem.align` -- Aligns the given reads with BWA mem.


The options for running BWA can be customized via the three Options classes:
    - :class:`~samwell.bwa_mem.AlgorithmOptions` -- Bwa mem algorithm options.
    - :class:`~samwell.bwa_mem.ScoringOptions` -- Bwa mem scoring options
    - :class:`~samwell.bwa_mem.InputOutputOptions` -- Bwa mem input and output options

An input read for alignment must be minimally transformed into a FASTQ-like record:

    - :class:`~samwell.bwa_mem.FastqRecord` -- Fastq record used as input to alignment.


Implementation
~~~~~~~~~~~~~~

Alignment of reads is performed asynchronously.

This is achieved by creating three sub-processes in :func:`~samwell.bwa_mem.align`:
    1. A process to consume the input iterable of reads and write them (a) to the stdin of
    BWA mem, and (b) to a queue of reads that are awaiting alignment results from BWA mem.
    2. A process to run BWA mem, where FASTQ records are written to the process' stdin,
    alignment results are returned to stdout, and any error/logging information from BWA mem is
    returned to stderr.
    3. A process to route the stderr of BWA mem to the given stderr handle (stderr_out).

Then :func:`~samwell.bwa_mem.align` method consumes the stdout of the BWA mem process and collates
that with the queue of reads that have been written/given to BWA mem (from process 1).  For each
input read, one or more alignments is expected to be returned by BWA mem.  The order in which
alignments of reads are returned by BWA mem is the same order as the order of reads given to BWA
mem.  The :func:`~samwell.bwa_mem.align` method then returns an iterable over the alignment
results.

Exceptions may occur in the thread to input FASTQ records to BWA mem, which are propagated
to the caller.  Furthermore, an exception is returned if the # of reads given to BWA mem is
not the same as the # of reads returned by BWA mem.

Some specific handling occurs around reading the BWA mem output with :py:mod:`~pysam`, since the
latter blocks waiting for at least some reads from BWA mem, which may not happen if there was an
issue in the various upstream processes (input to BWA mem or BWA mem itself).  This would have
caused a deadlock.


Examples
~~~~~~~~

Typically, we have :class:`~pysam.AlignedSegment` records obtained from reading from a SAM or BAM
file.  The first must be converted into :class:`~samwell.bwa_mem.FastqRecord` objects.

.. code-block:: python

    >>> from samwell.sam.bwa_mem import FastqRecord
    >>> from samwell.sam import reader
    >>> reads = reader("/path/to/sample.sam")
    >>> fastq_reads = map(lambda read: FastqRecord.build(read), reads)

Next, those :class:`~samwell.bwa_mem.FastqRecord`s can be aligned with the
:func:`~samwell.bwa_mem.align_mem` method.

.. code-block:: python

    >>> from samwell.sam.bwa_mem import align
    >>> results = map(lambda read: align(read), fastq_reads)

This returns an iterable over the alignment results.  An alignment result is a tuple
consisting of the original :class:`~samwell.bwa_mem.FastqRecord` and an iterator over the
alignments (see :class:`~pysam.AlignedSegment`).

.. code-block:: python

    >>> result = next(result)
    >>> fastq_read, alignments = result
    >>> str(fastq_read)
    @name
    GATTACA
    +
    HIJKLKM
    >>> len(alignments)
    2
    >>> alignment = str(next(alignments))
    >>> alignment.query_name
    name
    >>> type(alignment)
    <class 'pysam.libcalignedsegment.AlignedSegment'>
"""


import enum
import logging
import queue
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Any
from typing import Callable
from typing import ClassVar
from typing import Dict
from typing import Generic
from typing import Iterable
from typing import Iterator
from typing import List
from typing import Optional
from typing import Tuple
from typing import TypeVar
from typing import Union
from typing import cast

import attr
import pysam

import samwell.sam as sam
from samwell.itertools import PeekableIterator
from samwell.sam import SamFileType
from samwell.dnautils import reverse_complement

# The type for the source attribute for a :class:`samwell.bwa_mem.FastqRecord`
FastqRecordSourceType = TypeVar('FastqRecordSourceType')


@attr.s(frozen=True, auto_attribs=True)
class FastqRecord:
    """Fastq record used as input to alignment.

      Attributes:
          name: the name of the read
          bases: the read bases
          quals: the base qualities
          source: optionally the :class:`~pysam.AlignedSegment` from which this was built
          needs_alignment: True if the read needs alignment, False otherwise
          read_number: optionally the read number; should be set to 1 or 2 for paired end
              reads.
      """

    name: str = attr.ib()
    bases: str = attr.ib()
    quals: str = attr.ib()
    source: Optional[FastqRecordSourceType] = None  # type: ignore
    needs_alignment: bool = True
    read_number: Optional[int] = None

    _BASE_QUALITY_OFFSET: ClassVar[int] = 33

    @classmethod
    def build(cls,
              read: pysam.AlignedSegment,
              needs_alignment: bool = True,
              aligned_bases_only: bool = False,
              clip_three_prime: int = 0
              ) -> 'FastqRecord':
        """Builds a :class:`~samwell.bwa_mem.FastqRecord` from a :class:`~pysam.AlignedSegment`

        Args:
            read: the read to convert
            needs_alignment: True if the read should be aligned, False otherwise
            aligned_bases_only: only align the aligned bases (excludes soft-clipped bases)
            clip_three_prime: the number of bases to clip on the three-prime end of the read
                relative to the original direction of sequencing.  This will be applied after
                extracting the bases based on ``aligned_bases_only``.
        """
        # Get the bases and qualities
        if needs_alignment:
            if aligned_bases_only:
                bases = read.query_alignment_sequence
                quals = read.query_alignment_qualities
            else:
                bases = read.query_sequence
                quals = read.query_qualities

            # reverse complement if necessary
            if read.is_reverse:
                bases = reverse_complement(bases)
                quals = quals[::-1]

            if clip_three_prime > 0:
                index_from_end = -1 * clip_three_prime
                bases = bases[:index_from_end]
                quals = quals[:index_from_end]

            # convert to string
            quals = "".join([chr(q + FastqRecord._BASE_QUALITY_OFFSET) for q in quals])
        else:
            # If we're not going to align it, no need to muck with bases and quals
            bases = ""
            quals = ""

        # Get the read number
        if read.is_paired:
            read_number = 1 if read.is_read1 else 2
        else:
            read_number = None

        return FastqRecord(name=read.query_name,
                           bases=bases,
                           quals=quals,
                           source=read,
                           needs_alignment=needs_alignment,
                           read_number=read_number)

    def __hash__(self) -> int:
        """Returns a unique value for this record given the inputs.

        If source is defined and is a :class:`~pysam.AlignedSegment`, then the source's hash
        will be returned.  Otherwise, the hash of the concatenation of the name, bases, and
        qualities will be returned.
        """
        if self.source is not None and issubclass(self.source, pysam.AlignedSegment):
            return hash(self.source)
        else:
            return hash(self.str_with_read_number())

    def __str__(self) -> str:
        return f"@{self.name}\n{self.bases}\n+\n{self.quals}\n"

    def str_with_read_number(self) -> str:
        """Returns the record in FASTQ format, with the read number appended (colon delimited)."""
        name = self.name + ":" + (str(self.read_number) if self.read_number is not None else "0")
        return f"@{name}\n{self.bases}\n+\n{self.quals}\n"


class _CommandLineOptionGroup:
    """Base class for groups of bwa options using @attr.s.

    It is assumed that every attribute has the 'flag' key specified in its metadata field.  Use the
    :func:`~samwell.bwa_mem._flag` method to add additional flag attributes.
    """

    def args(self) -> List[str]:
        """Build the list of command line arguments from the defined options."""
        _args = []
        flag_to_attribute_name: Dict[str, str] = {}
        # go through each attribute
        for attribute in attr.fields(type(self)):
            # get the value for the flag
            value = getattr(self, attribute.name)
            if isinstance(value, enum.Enum):
                value = value.value
            else:
                # check if it iterable, and if so, join them with commas
                try:
                    value = ",".join(iter(value))
                except TypeError:
                    pass

            # if it is set, add it to args
            if value is not None:
                # assume that they have metadata, with the "flag" specified.  Get the flag to use
                flag = attribute.metadata['flag']
                if flag in flag_to_attribute_name:
                    cur_name = attribute.name
                    other_name = flag_to_attribute_name[flag]
                    raise ValueError(
                        f"Flag '{flag}' found in attributes {cur_name} and {other_name}")
                flag_to_attribute_name[flag] = attribute.name
                if attribute.type in (bool, Optional[bool]):
                    if value is True:
                        _args.append(flag)
                else:
                    _args.extend([flag, str(value)])
        return _args


# Alias for the alignment result
AlignmentResult = Tuple[FastqRecord, List[pysam.AlignedSegment]]


@attr.s(frozen=True)
class AlgorithmOptions(_CommandLineOptionGroup):
    """Bwa mem algorithm options

    Attributes:
        threads: number of threads
        min_seed_len: minimum seed length
        band_width: band width for banded alignment
        off_diagonal_dropoff: off-diagonal X-dropoff
        internal_seeds_length_factor: look for internal seeds inside a seed longer than
            min_seed_len * internal_seeds_length_factor
        max_third_seed_occurrence: seed occurrence for the 3rd round seeding
        max_seed_occurrence: skip seeds with more than INT occurrences
        drop_ratio: drop chains shorter than this fraction of the longest overlapping chain
        min_chain_weight: discard a chain if seeded bases shorter than this value
        max_mate_rescue_rounds: perform at most INT rounds of mate rescues for each read
        skip_mate_rescue: skip mate rescue
        skip_pairing: skip pairing; mate rescue performed unless skip_mate_rescue also in use
    """

    threads: Optional[int] = attr.ib(default=None, metadata={'flag': '-t'})
    min_seed_len: Optional[int] = attr.ib(default=None, metadata={'flag': '-k'})
    band_width: Optional[int] = attr.ib(default=None, metadata={'flag': '-w'})
    off_diagonal_dropoff: Optional[int] = attr.ib(default=None, metadata={'flag': '-d'})
    internal_seeds_length_factor: Optional[float] = attr.ib(default=None, metadata={'flag': '-r'})
    max_third_seed_occurrence: Optional[int] = attr.ib(default=None, metadata={'flag': '-y'})
    max_seed_occurrence: Optional[int] = attr.ib(default=None, metadata={'flag': '-c'})
    drop_ratio: Optional[float] = attr.ib(default=None, metadata={'flag': '-D'})
    min_chain_weight: Optional[int] = attr.ib(default=None, metadata={'flag': '-W'})
    max_mate_rescue_rounds: Optional[int] = attr.ib(default=None, metadata={'flag': '-m'})
    skip_mate_rescue: Optional[bool] = attr.ib(default=None, metadata={'flag': '-S'})
    skip_pairing: Optional[bool] = attr.ib(default=None, metadata={'flag': '-P'})


@enum.unique
class ReadType(enum.Enum):
    """The read type for BWA mem."""

    PacBio = "pacbio"
    OxfordNano2D = "ont2d"
    IntraSpecies = "intractg"


@attr.s(frozen=True)
class ScoringOptions(_CommandLineOptionGroup):

    """Bwa mem scoring options

    Attributes:
        match_score: the score for a sequence match, which scales options -TdBOELU unless
            overridden
        mismatch_score: penalty for a mismatch
        gap_open: gap open penalties for deletions and insertions (single value to use the same
            for both)
        gap_extend: gap extension penalty; a gap of size k cost '{-O} + {-E}*k' (single value to
            use the same for both)
        clipping_penalty: penalty for 5'- and 3'-end clipping
        unpaired_penalty: penalty for an unpaired read pair
        read_type: read type. Setting -x changes multiple parameters unless overriden:
                     pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)
                     ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)
                     intractg: -B9 -O16 -L5  (intra-species contigs to ref)
    """

    match_score: Optional[int] = attr.ib(default=None, metadata={'flag': '-A'})
    mismatch_score: Optional[int] = attr.ib(default=None, metadata={'flag': '-B'})
    gap_open: Optional[Union[int, Tuple[int, int]]] = attr.ib(default=None,
                                                              metadata={'flag': '-O'})
    gap_extend: Optional[Union[int, Tuple[int, int]]] = attr.ib(default=None,
                                                                metadata={'flag': '-E'})
    clipping_penalty: Optional[Union[int, Tuple[int, int]]] = attr.ib(default=None,
                                                                      metadata={'flag': '-L'})
    unpaired_penalty: Optional[int] = attr.ib(default=None, metadata={'flag': '-U'})
    read_type: Optional[ReadType] = attr.ib(default=None, metadata={'flag': '-x'})


# The type for BWA mem's insert size parameter (-I) option
InsertSizeParamsType = Union[
    float,
    Tuple[float, float],
    Tuple[float, float, int],
    Tuple[float, float, int, int]]


@attr.s(frozen=True)
class InputOutputOptions(_CommandLineOptionGroup):
    """Bwa mem input and output options

    Attributes:
        interleaved_pairs: read pairs are consecutive (r1 then r2), otherwise fragment reads
        read_group: read group header line such as '@RG\tID:foo\tSM:bar
        header_insert: insert STR to header if it starts with @; or insert lines in FILE
        alts_as_primary: treat ALT contigs as part of the primary assembly (i.e. ignore
            <idxbase>.alt file)
        verbosity: verbose level: 1=error, 2=warning, 3=message, 4+=debuggin
        min_alignment_score: minimum score to output
        max_hits_within_max_score: if there are <INT hits with score >80% of the max score, output
            all in XA
        all_alignments: output all alignments for SE or unpaired PE
        append_fastq_comment: append FASTA/FASTQ comment to SAM output
        add_fasta_header_to_xr: output the reference FASTA header in the XR tag
        softclip_supplementary: use soft clipping for supplementary alignments
        split_hits_are_secondary: mark shorter split hits as secondary
        insert_size_params: specify the mean, standard deviation (10% of the mean if absent), max
            (4 sigma from the mean if absent) and min of the insert size distribution.  FR
            orientation only.
        bases_per_batch: how many bases of sequence data bwa should read from the input
            before triggering a batch of alignments.
    """

    interleaved_pairs: Optional[bool] = attr.ib(default=None, metadata={'flag': '-p'})
    read_group: Optional[str] = attr.ib(default=None, metadata={'flag': '-R'})
    header_insert: Optional[Union[str, Path]] = attr.ib(default=None, metadata={'flag': '-H'})
    alts_as_primary: Optional[bool] = attr.ib(default=None, metadata={'flag': '-j'})
    verbosity: Optional[int] = attr.ib(default=None, metadata={'flag': '-v'})
    min_alignment_score: Optional[int] = attr.ib(default=None, metadata={'flag': '-T'})
    max_hits_within_max_score: Optional[Union[int, Tuple[int, int]]] = \
        attr.ib(default=None, metadata={'flag': '-h'})
    all_alignments: Optional[bool] = attr.ib(default=None, metadata={'flag': '-a'})
    append_fastq_comment: Optional[bool] = attr.ib(default=None, metadata={'flag': '-C'})
    add_fasta_header_to_xr: Optional[bool] = attr.ib(default=None, metadata={'flag': '-V'})
    softclip_supplementary: Optional[bool] = attr.ib(default=None, metadata={'flag': '-Y'})
    split_hits_are_secondary: Optional[bool] = attr.ib(default=None, metadata={'flag': '-M'})
    insert_size_params: Optional[InsertSizeParamsType] = attr.ib(default=None,
                                                                 metadata={'flag': '-I'})
    bases_per_batch: Optional[int] = attr.ib(default=115000, metadata={'flag': '-K'})


# The type for the source items in :class:`samwell.bwa_mem._SourceToSinkThread`
SourceToSinkThreadType = TypeVar('SourceToSinkThreadType')


class _SourceToSinkThread(threading.Thread, Generic[SourceToSinkThreadType]):
    """A thread that consumes elements from the source and adds them to the sink

    Attributes:
        num_added: the number of elements from the source added to the sink
    """

    def __init__(self,
                 source: Iterator[SourceToSinkThreadType],
                 sink_add_func: Callable[[SourceToSinkThreadType], None],
                 sink_close_func: Optional[Callable[..., None]] = None) -> None:
        """Creates a new thread for consuming the source and adding to the sink.

        Args:
            source: the source iterator from which to consume
            sink_add_func: the method to use to add an element to the sink
            sink_close_func: the method to call when all elements have been added to the sink
        """
        super().__init__(daemon=True)
        self.num_added: int = 0
        self._source = source
        self._sink_add_method = sink_add_func
        self._sink_close_method = sink_close_func
        self.exception: Optional[Exception] = None
        self.done = False

    def run(self) -> None:
        """Runs the source to sink transfer"""
        try:
            for item in self._source:
                self._sink_add_method(item)
                self.num_added += 1
        except Exception as e:
            self.exception = e
        finally:
            if self._sink_close_method is not None:
                self._sink_close_method()
            self.done = True


def _same_read(read: FastqRecord, alignment: pysam.AlignedSegment) -> bool:
    """True if an alignment for the given read, False otherwise.

    For the alignment to be considered as an alignment for this read, the read name and read number
    must match.  The read number is appended to the alignment's query name.
    """
    if alignment.is_paired:
        assert read.read_number is not None, f"Paired alignment but the read has no read #: {read}"
        if read.name != alignment.query_name:
            return False
        else:
            alignment_read_number = 1 if alignment.is_read1 else 2
            return read.read_number == alignment_read_number
    else:
        alignment_name, alignment_read_number = alignment.query_name.rsplit(':', 1)
        if read.name != alignment_name:
            return False
        elif read.read_number is None:
            return int(alignment_read_number) == 0
        else:
            return read.read_number == int(alignment_read_number)


def _collate_alignments(reads_queue: queue.Queue,
                        alignments_reader: pysam.AlignmentFile,
                        suppress_secondaries: bool = False) -> Iterable[AlignmentResult]:
    """Collates the alignments for each read in the given queue.

    Alignments for reads in the alignments reader are in the same order as the reads in the reads
    queue.  This allows traversal of both once and in a step-wise fashion.  In fact, there exist
    alignments for reads in the read queue that need alignment (see the corresponding property).

    This method may block waiting for (1) reads to be written for BWA to consume and thus added
    to the queue that will be consumed by this method, or (2) the alignments returned by BWA for a
    given read.  The former (1) will not block indefinitely since at some point the BWA mem
    input process will write the sentinel value and reads_queue.get will return None.  The
    latter will not block indefinitely since the stdout will be closed when the BWA mem process is
    terminated.

    Args:
        reads_queue: the queue of reads
        alignments_reader: the reader of sam records
        suppress_secondaries: true to discard all secondary alignments, false otherwise

    Returns:
        An iterable over the alignment results.  An alignment result is a tuple consisting of the
        original :class:`~samwell.bwa_mem.FastqRecord` and an iterator over the alignments (see
        :class:`~pysam.AlignedSegment`).
    """
    alignments_iter: PeekableIterator = PeekableIterator(alignments_reader)
    reads_iterator = cast(Iterator[FastqRecord], iter(reads_queue.get, None))
    for read in reads_iterator:
        results: List[pysam.AlignedSegment] = []

        if read.needs_alignment:
            result = alignments_iter.peek() if alignments_iter.can_peek() else None
            while result is not None and _same_read(read, result):
                next(alignments_iter)  # consume the current record
                if not suppress_secondaries or not result.is_secondary:
                    # Update the query name since we may have originally appended the read number
                    result.query_name = read.name
                    results.append(result)
                result = alignments_iter.peek() if alignments_iter.can_peek() else None

        yield (read, results)
    assert not alignments_iter.can_peek(), 'Alignments exist but no more reads in the queue'


def _build_command_line(idxbase: Path,
                        executable_path: Path = Path('bwa'),
                        algo_opts: Optional[AlgorithmOptions] = None,
                        scoring_opts: Optional[ScoringOptions] = None,
                        io_opts: Optional[InputOutputOptions] = None) -> List[str]:
    """Builds the command line for bwa mem.

    Args:
        idxbase: the path prefix for all the BWA-specific index files
        executable_path: the path to the BWA executable
        algo_opts: the algorithm options
        scoring_opts: the scoring options
        io_opts: the input and output options
    """
    # Start with the path to BWA, then the mem command
    cmd: List[Any] = [executable_path, 'mem']
    # Add any options
    for opts in [algo_opts, scoring_opts, io_opts]:
        if opts is not None:
            cmd.extend(opts.args())
    # Now the reference genome index basename
    cmd.append(idxbase)
    # Now set the input to be from standard input
    cmd.append('/dev/stdin')
    # Convert all args to strings
    args = [str(arg) for arg in cmd]
    return args


def _build_bwa_input_process(reads: Iterable[FastqRecord],
                             to_bwa_handle: Any,
                             to_output_queue: queue.Queue,
                             interleaved_pairs: Optional[bool] = None
                             ) -> _SourceToSinkThread:
    """Builds and starts a process to write the given FASTQ records for BWA mem and the given queue

    Args:
        reads: the reads to input to the BWA mem subprocess' stdin
        to_bwa_handle: the IO handle to which to write the FASTQ records for BWA mem
        to_output_queue: the queue to also write to after writing a FASTQ record for BWA mem;
            this queue is mainly used for collating FASTQ reads and SAM alignments.
        interleaved_pairs: read pairs are consecutive (r1 then r2), otherwise unpaired reads
    """
    last_read_name_and_number: Optional[Tuple[str, Optional[int]]] = None

    def sink_add_method(read: FastqRecord) -> None:
        """Writes a FASTQ record to the BWA mem input as well as the results collation"""
        nonlocal last_read_name_and_number

        if read.needs_alignment:
            if last_read_name_and_number is not None:
                name, read_number = last_read_name_and_number
                assert name != read.name or read_number != read.read_number, \
                    'Consecutive reads have the same name and read number:' + \
                    f'\n\t\tname: {name}\n\t\tread number: {read_number}' + \
                    f'\n\t\tsource: {read.source}'
            last_read_name_and_number = (read.name, read.read_number)
            if interleaved_pairs is True:
                to_bwa_handle.write(str(read))
            else:
                # IMPORTANT: the read name has the read number appended to disambiguate ends of a
                # pair
                to_bwa_handle.write(read.str_with_read_number())
        to_output_queue.put(read)

    def sink_close_method() -> None:
        """Close the BWA mem input handle and the output queue for results collation"""
        to_bwa_handle.close()
        to_output_queue.put(None)  # add the sentinel value

    bwa_input_process = _SourceToSinkThread(source=iter(reads),
                                            sink_add_func=sink_add_method,
                                            sink_close_func=sink_close_method)
    bwa_input_process.start()

    return bwa_input_process


def align(reads: Iterable[FastqRecord],
          idxbase: Path,
          executable_path: Path = Path('bwa'),
          algo_opts: Optional[AlgorithmOptions] = None,
          scoring_opts: Optional[ScoringOptions] = None,
          io_opts: Optional[InputOutputOptions] = None,
          suppress_secondaries: bool = False,
          stderr_out: Any = sys.stderr
          ) -> Iterable[AlignmentResult]:
    """Aligns the given reads with BWA mem.

    See :py:mod:`~samwell.bwa_mem` for a detailed explanation for the implementation approach.

    Args:
        reads: the reads to align
        idxbase: the path prefix for all the BWA-specific index files
        executable_path: the path to the BWA executable
        algo_opts: the algorithm options
        scoring_opts: the scoring options
        io_opts: the input and output options
        suppress_secondaries: true to discard all secondary alignments, false otherwise

    Returns:
        An iterable over the alignment results.  An alignment result is a tuple consisting of the
        original :class:`~samwell.bwa_mem.FastqRecord` and an iterator over the alignments (see
        :class:`~pysam.AlignedSegment`)
    """

    # Build the command line used to run BWA MEM
    command_line = _build_command_line(idxbase=idxbase,
                                       executable_path=executable_path,
                                       algo_opts=algo_opts,
                                       scoring_opts=scoring_opts,
                                       io_opts=io_opts)

    # Create a sub-process in which to run BWA mem.  This process will read FASTQ records from
    # stdin, write SAM records to stdout, and write any error/logging information to stderr.
    bwa_mem_process = subprocess.Popen(args=command_line,
                                       stdin=subprocess.PIPE,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       universal_newlines=True)

    # Create a sub-process in which we read the stderr of the BWA mem subprocess and write it to
    # the given stderr_out handle.
    bwa_mem_stderr_process = _SourceToSinkThread(source=iter(bwa_mem_process.stderr),
                                                 sink_add_func=stderr_out.write,
                                                 sink_close_func=None)
    bwa_mem_stderr_process.start()

    # Create a queue of FASTQ records that the sub-process who will write to BWA mem's stdin
    # will also write.  This is so we can collate/join the input FASTQ records with the output SAM
    # (or alignment) records.  A sentinel value (None) will be written to indicate no more reads
    # will be placed in the queue.
    reads_queue: queue.Queue = queue.Queue()

    # Create a sub-process to consume the input FASTQ records and write them to BWA mem's stdin. We
    # write in a separate thread to avoid any deadlock with waiting for output from BWA mem's
    # stdout.  This can happen in a synchronous implementation where BWA mem is buffering reads and
    # we are waiting for some results from BWA mem's stdout, but really BWA mem is waiting for
    # either more reads from stdin or for stdin to be closed.
    interleaved_pairs = io_opts.interleaved_pairs if io_opts is not None else None
    bwa_input_process = _build_bwa_input_process(reads=reads,
                                                 to_bwa_handle=bwa_mem_process.stdin,
                                                 to_output_queue=reads_queue,
                                                 interleaved_pairs=interleaved_pairs)

    # Go through the output
    num_aligned = 0
    try:
        # Wait for some reads to be written.  pysam will block opening the input file until some
        # data is available, or the stream is closed.  If no data is added, don't even try opening
        # the stream.
        while bwa_input_process.num_added == 0 and not bwa_input_process.done:
            # the input process is still running but no reads have been added
            time.sleep(.1)
        if bwa_input_process.num_added == 0 and bwa_input_process.done:
            # the input process is done (error or success) and no reads have been added, so skip
            # opening pysam
            raise StopIteration
        # Read through the output of BWA mem, and collate that with the queue of reads given to
        # BWA mem
        with sam.reader(path=bwa_mem_process.stdout, file_type=SamFileType.SAM) as reader:
            alignment_results = _collate_alignments(reads_queue=reads_queue,
                                                    alignments_reader=reader,
                                                    suppress_secondaries=suppress_secondaries)
            # A simple loop with its only purpose to count the number of alignment results
            for result in alignment_results:
                num_aligned += 1
                yield result
    finally:
        # Close the stdin of the BWA mem process.  This should signal BWA mem to shut down, and
        # for the input thread to stop.
        bwa_mem_process.stdin.close()

        # Join the input thread as now stdin of the BWA mem process is closed.
        bwa_input_process.join(timeout=1.0)

        # Check if the inputting reads to BWA had an exception
        if bwa_input_process.exception is not None:
            raise bwa_input_process.exception
        elif bwa_input_process.is_alive():
            raise RuntimeError("BWA process encountered no errors but did not terminate.")

        # Check that the number of reads given to BWA mem was the same # returned by BWA mem
        num_left = bwa_input_process.num_added - num_aligned
        if num_left != 0:
            raise ValueError(f"Still had {num_left:,d} remaining reads from BWA")

        # Shut down the BWA mem process.  If it fails to shutdown, log a warning and continue on
        try:
            bwa_mem_process.wait(timeout=5.0)
        except subprocess.TimeoutExpired as ex:
            logger = logging.getLogger(__name__)
            logger.warning("Could not shutdown BWA, ignoring error: %s", str(ex))

        # Shut down the stderr thread
        bwa_mem_stderr_process.join(timeout=1.0)
