import unittest

import yaml

# Adjust the import according to your module structure
from frmatcher.fastq_file_name_checker import FastqFileNameChecker


class TestFastqFileNameChecker(unittest.TestCase):
    def test_categorization(self):
        # Define the YAML patterns directly in the test
        patterns = yaml.safe_load(
            r"""
        patterns:
          r1:
            - "_1"
            - "_R1"
          r2:
            - "_2"
            - "_R2"
          ignore:
            - "^i_"
            - "^I_"
            - "_i\\d+"
            - "_I\\d+"
        """
        )

        filenames = [
            "sample_1_L001.fastq.gz",  # R1
            "sample_R1_L001.fastq.gz",  # R1
            "sample_2_L001.fastq.gz",  # R2
            "sample_R2_L001.fastq.gz",  # R2
            "sample_i1_L001.fastq.gz",  # Ignored
            "sample_I2_L001.fastq.gz",  # Ignored
            "i_sample_1_L001.fastq.gz",  # Ignored
            "I_sample_2_L001.fastq.gz",  # Ignored
            "sample_A_L001.fastq.gz",  # Ignored (no matching pattern)
        ]

        # Initialize checker with in-memory patterns
        checker = FastqFileNameChecker(filenames, config_path=None)
        checker.patterns = patterns["patterns"]  # Inject the test patterns directly
        categorized_files = checker.categorize_fastq_files()

        # Assert correct categorization
        expected_r1 = {"sample_1_L001.fastq.gz", "sample_R1_L001.fastq.gz"}
        expected_r2 = {"sample_2_L001.fastq.gz", "sample_R2_L001.fastq.gz"}
        expected_ignored = {
            "sample_i1_L001.fastq.gz",
            "sample_I2_L001.fastq.gz",
            "i_sample_1_L001.fastq.gz",
            "I_sample_2_L001.fastq.gz",
            "sample_A_L001.fastq.gz",
        }

        with self.subTest("R1 Categorization"):
            self.assertEqual(set(categorized_files["R1"]), expected_r1)

        with self.subTest("R2 Categorization"):
            self.assertEqual(set(categorized_files["R2"]), expected_r2)

        with self.subTest("Ignored Categorization"):
            self.assertEqual(set(categorized_files["ignored"]), expected_ignored)


if __name__ == "__main__":
    unittest.main()
