![Python package](https://github.com/myriad-opensource/samwell/workflows/Python%20package/badge.svg)

# Samwell: a python package for using genomic files... well

Samwell provides elegant utilities for managing biological data.

## Quickstart

First install samwell:

```
pip install samwell
```

### Reading/Writing BAMs with automatic inference of filetype

Samwell provides easy utilities for reading/writing BAMs:

```python
from samwell import sam
with sam.reader("myfile.bam") as in_bam:
    with sam.writer("my-output-file.bam", header=in_bam.header) as out_bam:
        for read in in_bam:
            if read.is_paired:
                out_bam.write(read)
```


### Realigning fastqs with bwa

You can use `samwell` to easily realign fastq records as necessary

```python
from pathlib import Path
from samwell import sam
from samwell.sam import bwa_mem
from samwell.sam import clipping
from samwell.sam.bwa_mem import FastqRecord
with sam.reader("myfile.bam") as in_bam:
    with sam.writer("outfile.bam", header=in_bam.header) as out_bam:
         fastq_gen = iter(FastqRecord.build(read) for read in in_bam)
         for read in bwa_mem.align(fastq_gen, Path("genome.fasta")):
             out_bam.write(read)
```

See `samwell.bwa_mem` module for more detail.


## Developing with samwell

Samwell uses [`poetry`](https://github.com/python-poetry/poetry#installation) for dependency managment.

Please install `poetry` using the instructions in the above link.
Then simply execute:

```bash
poetry install
```

## Checking the Build

### Linting 

```bash
poetry run flake8 --config=flake8.cfg samwell
```

### Type Checking

```bash
poetry run mypy -p samwell --config=mypy.ini
```

### Unit Tests

```bash
poetry run python -m pytest --cov=samwell --cov-branch
```
