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
with sam.writer("my-output-file.bam") as out_bam:
    with sam.reader("myfile.bam") as in_bam:
        for read in fp:
            if read.is_paired:
                out_bam.write(read)
```


### Realigning fastqs with bwa

You can use `samwell` to easily realign fastq records as necessary

```python
from samwell import sam
from samwell import bwa_mem
with sam.reader("myfile.bam") as in_bam:
    with sam.writer("outfile.bam") as out_bam:
        fastq_gen = (FastqRecord.build(read) for read in in_bam)
        gen = clip_specific_reads(fastq_gen)
        out_bam.write(bwa_mem.align(gen))
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
