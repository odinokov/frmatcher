# FRMatcher

**FRMatcher** categorizes a list of presumably FASTQ files into `R1` (forward reads) and `R2` (reverse reads) pairs using customizable pattern matching.

## Installation

Clone the repository:
   ```bash
   git clone https://github.com/odinokov/frmatcher.git
   cd frmatcher
   ```

Activate the virtual environment:
   ```bash
   poetry shell
   ```

Build the package:
   ```bash
   poetry build
   ```

Install the package locally:
   ```bash
   poetry install
   ```

## Usage

```python
from frmatcher import FastqFileNameChecker

filenames = [
    "sample_1_L001.fastq.gz",
    "sample_2_L001.fastq.gz",
    "sample_1_L002.fastq.gz",
    "sample_2_L002.fastq.gz",
]

checker = FastqFileNameChecker(filenames,
                              length_check=False,
                              verbose=False)

# checker = FastqFileNameChecker(filenames,
#                                length_check=True,
#                                verbose=True,
#                                config_path=None)
# checker.patterns = {
#     'r1': ["_1", "_R1"],
#     'r2': ["_2", "_R2"],
#     'ignore': ["^i_", "^I_", "_i\\d+", "_I\\d+"]
# }

categorized_files = checker.categorize_fastq_files()

print(categorized_files)

# {'R1': ['sample_1_L001.fastq.gz', 'sample_1_L002.fastq.gz'], 
# 'R2': ['sample_2_L001.fastq.gz', 'sample_2_L002.fastq.gz'], 
# 'ignored': []}

```

## License

MIT License. See the [LICENSE](LICENSE) file for details.
