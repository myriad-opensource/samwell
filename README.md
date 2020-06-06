[![Language][language-badge]][language-link]
[![Code Style][code-style-badge]][code-style-link]
[![Type Checked][type-checking-badge]][type-checking-link]
[![PEP8][pep-8-badge]][pep-8-link]
[![License][license-badge]][license-link]

---

[![Python package][python-package-badge]][python-package-link]
[![PyPI version][pypi-badge]][pypi-link]
[![PyPI download total][pypi-downloads-badge]][pypi-downloads-link]

---

[language-badge]:       http://img.shields.io/badge/language-python-brightgreen.svg
[language-link]:        http://www.python.org/
[code-style-badge]:     https://img.shields.io/badge/code%20style-black-000000.svg
[code-style-link]:      https://black.readthedocs.io/en/stable/ 
[type-checking-badge]:  http://www.mypy-lang.org/static/mypy_badge.svg
[type-checking-link]:   http://mypy-lang.org/
[pep-8-badge]:          https://img.shields.io/badge/code%20style-pep8-brightgreen.svg
[pep-8-link]:           https://www.python.org/dev/peps/pep-0008/
[license-badge]:        http://img.shields.io/badge/license-MIT-blue.svg
[license-link]:         https://github.com/myriad-opensource/samwell/blob/master/LICENSE
[python-package-badge]: https://github.com/myriad-opensource/samwell/workflows/Python%20package/badge.svg
[python-package-link]:  https://github.com/myriad-opensource/samwell/actions?query=workflow%3A%22Python+package%22
[pypi-badge]:           https://badge.fury.io/py/samwell.svg
[pypi-link]:            https://pypi.python.org/pypi/samwell
[pypi-downloads-badge]: https://img.shields.io/pypi/dm/samwell
[pypi-downloads-link]:  https://pypi.python.org/pypi/samwell

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
