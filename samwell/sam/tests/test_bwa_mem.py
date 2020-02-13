import distutils.spawn
import subprocess
from pathlib import Path
from tempfile import NamedTemporaryFile as NamedTemp
from typing import List
from typing import Optional
from typing import Tuple

import attr
import pytest
from py._path.local import LocalPath as TmpDir
from pysam import AlignedSegment

from samwell.sam.bwa_mem import FastqRecord
from samwell.sam.bwa_mem import InputOutputOptions
from samwell.sam.bwa_mem import align


BwaExecutable: Optional[str] = distutils.spawn.find_executable("bwa")


@pytest.fixture
def ref_fasta(tmpdir: TmpDir) -> Path:
    with NamedTemp(suffix=".fasta", dir=tmpdir, mode='w', delete=False) as fp:
        filename = Path(fp.name).name
    ref_fasta = tmpdir / filename

    with ref_fasta.open('w') as fh:
        fh.write(">1\n")
        fh.write("CCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACCCTAACCCTAA\n")
    subprocess.check_call(args=["bwa", "index", str(fp.name)])
    return ref_fasta


@pytest.fixture
def fastq_record() -> FastqRecord:
    read_bases = "CCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACCCTAACCCTAA"
    return FastqRecord(
        name="some_name",
        bases=read_bases,
        quals="".join(["I" for _ in read_bases])
    )


def _assert_alignment_for_fastq_record(read: FastqRecord,
                                       results: Tuple[FastqRecord, List[AlignedSegment]]) -> None:

    fastq, alignments = results
    assert len(alignments) == 1
    alignment = alignments[0]

    assert alignment.query_name == read.name
    assert alignment.query_sequence == read.bases
    assert "".join([chr(q + 33) for q in alignment.query_qualities]) == read.quals
    assert alignment.query_name == read.name
    assert alignment.reference_name == "1"
    assert alignment.reference_start == 0
    assert alignment.cigarstring == "60M"
    if read.read_number is not None:
        assert read.read_number == 1 or read.read_number == 2
        assert alignment.is_paired
        assert read.read_number == (1 if alignment.is_read1 else 2)


@pytest.mark.skipif(BwaExecutable is None, reason="requires bwa 0.7.17")
def test_single_alignment(fastq_record: FastqRecord, ref_fasta: Path) -> None:
    # run BWA
    results = list(align(reads=[fastq_record], idxbase=ref_fasta))

    # Check the returned alignments
    assert len(results) == 1
    _assert_alignment_for_fastq_record(read=fastq_record, results=results[0])


@pytest.mark.skipif(BwaExecutable is None, reason="requires bwa 0.7.17")
def test_fails_consecutive_reads_with_the_same_name_and_number(fastq_record: FastqRecord,
                                                               ref_fasta: Path) -> None:
    # run BWA
    with pytest.raises(Exception, match="Consecutive reads"):
        list(align(reads=[fastq_record, fastq_record], idxbase=ref_fasta))


@pytest.mark.skipif(BwaExecutable is None, reason="requires bwa 0.7.17")
def test_paired_end_reads(fastq_record: FastqRecord, ref_fasta: Path) -> None:
    # run BWA
    r1 = attr.evolve(fastq_record, read_number=1)
    r2 = attr.evolve(r1, read_number=2)
    io_opts = InputOutputOptions(interleaved_pairs=True)
    results = list(align(reads=[r1, r2], idxbase=ref_fasta, io_opts=io_opts))
    _assert_alignment_for_fastq_record(read=r1, results=results[0])
    _assert_alignment_for_fastq_record(read=r2, results=results[1])


@pytest.mark.skipif(BwaExecutable is None, reason="requires bwa 0.7.17")
def test_needs_alignment(fastq_record: FastqRecord, ref_fasta: Path) -> None:
    # run BWA
    rec1 = fastq_record
    rec2 = FastqRecord(
        name="needs_alignment=False",
        bases=rec1.bases,
        quals=rec1.quals,
        needs_alignment=False
    )

    results = list(align(reads=[rec1, rec2], idxbase=ref_fasta))

    # Check the returned alignments
    assert len(results) == 2
    for result in results:
        fastq, alignments = result
        if fastq.needs_alignment:
            assert fastq == rec1
            _assert_alignment_for_fastq_record(read=fastq_record, results=result)
        else:
            assert len(alignments) == 0
            assert fastq.name == "needs_alignment=False"
            assert fastq == rec2


@pytest.mark.skipif(BwaExecutable is None, reason="requires bwa 0.7.17")
def test_no_alignment(ref_fasta: Path) -> None:
    fastq_record = FastqRecord(
        name="unmapped",
        bases="A" * 60,
        quals="I" * 60,
        needs_alignment=False
    )

    # run BWA
    results = list(align(reads=[fastq_record], idxbase=ref_fasta))

    # Check the returned alignments
    assert len(results) == 1
    fastq, alignments = results[0]
    assert len(alignments) == 0
    assert fastq_record == fastq
