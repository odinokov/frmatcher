# frmatcher

Scans a directory for FASTQ files and outputs a two-column TSV of paired R1/R2 paths.

Handles common naming conventions (`_R1_`, `_R2_`, `_1.`, `_2.`, `.R1.`, `.R2.`) and compound extensions (`.fastq.gz`, `.fq.bz2`, etc.). Index reads (`I1`/`I2`) are skipped.

## Installation

```bash
pip install git+https://github.com/odinokov/frmatcher.git
```

## CLI usage

```bash
frmatcher -i /path/to/fastq_dir -o pairs.tsv
```

```
options:
  -i, --in-dir    Directory containing FASTQ files (required)
  -o, --out-tsv   Output TSV file path (required)
  -r, --recursive Search subdirectories
  --strict        Exit non-zero if any non-index FASTQ cannot be paired or matched
  --log-level     TRACE|DEBUG|INFO|WARNING|ERROR|CRITICAL (default: INFO)
  --version       Show version and exit
```

### Output format

Each line is a paired R1/R2:

```
/data/sample1_R1_001.fastq.gz	/data/sample1_R2_001.fastq.gz
/data/sample2_R1_001.fastq.gz	/data/sample2_R2_001.fastq.gz
```

Output paths preserve the input directory style: absolute `--in-dir` values produce absolute paths, and relative `--in-dir` values produce relative paths.

## Python API

```python
from frmatcher import group_fastq_pairs, write_fastq_pairs, FASTQ_EXTENSIONS
from pathlib import Path
from loguru import logger

# Enable logging when using as a library
logger.enable("frmatcher")

# High-level: scan and write pairs in one call
write_fastq_pairs("/path/to/fastqs", "pairs.tsv", recursive=True)

# Low-level: inspect pairing results programmatically
stats = group_fastq_pairs(Path("/path/to/fastqs"), recursive=True)
print(f"Found {stats.n_fastq} FASTQ files, {len(stats.groups)} paired groups")
for key, reads in stats.groups.items():
    print(key, reads["1"], reads["2"])

# Custom extensions
from frmatcher import iter_fastq_files
custom_ext = FASTQ_EXTENSIONS | frozenset({".fq.zst"})
for path in iter_fastq_files(Path("/data"), extensions=custom_ext):
    print(path)
```

## Requirements

- Python >= 3.10
- [loguru](https://github.com/Delgan/loguru)
