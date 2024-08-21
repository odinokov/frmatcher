import re
from typing import Dict, List

import pkg_resources
import yaml
from loguru import logger


class FastqFileNameChecker:
    def __init__(
        self,
        filenames: List[str],
        config_path: str = None,
        length_check: bool = False,
        verbose: bool = False,
    ):
        """
        Initialize the FastqFileNameChecker with a list of filenames.

        Args:
            filenames (List[str]): List of filenames to categorize.
            config_path (str): Path to the YAML configuration file. Default loads from package if None.
            length_check (bool): Whether to check if all filenames have the same length. Default is False.
            verbose (bool): Whether to enable detailed logging. Default is False.

        Raises:
            ValueError: If filenames have different lengths and length_check is True.
        """
        self.filenames = filenames
        self.verbose = verbose

        # Set logging level based on verbosity
        logger.remove()  # Remove the default handler
        if self.verbose:
            logger.add(lambda msg: print(msg, end=""), level="DEBUG", colorize=True)
        else:
            logger.add(lambda msg: print(msg, end=""), level="ERROR", colorize=True)

        # Load patterns from the package if no config_path is provided
        if config_path is None:
            config_path = pkg_resources.resource_filename(
                __name__, "config/patterns.yaml"
            )
        self.patterns = self.load_patterns(config_path)

        if length_check:
            self._check_filename_lengths()

    def load_patterns(self, config_path: str) -> Dict[str, List[str]]:
        """
        Load patterns from a YAML configuration file.

        Args:
            config_path (str): Path to the YAML configuration file.

        Returns:
            Dict[str, List[str]]: A dictionary with R1, R2, and ignore patterns.
        """
        with open(config_path, "r") as file:
            config = yaml.safe_load(file)
        logger.debug(f"Loaded patterns from {config_path}")
        return config["patterns"]

    def _check_filename_lengths(self) -> None:
        """
        Checks if all filenames have the same length.

        Raises:
            ValueError: If filenames do not have the same length.
        """
        lengths = list(map(len, self.filenames))
        if len(set(lengths)) > 1:
            logger.error(
                "Filenames do not all have the same length. Please ensure all filenames are consistent."
            )
            raise ValueError(
                "Filenames do not all have the same length. Please ensure all filenames are consistent."
            )
        logger.info("All filenames have the same length.")

    def categorize_fastq_files(self) -> Dict[str, List[str]]:
        """
        Categorizes FASTQ files into R1, R2, or ignored based on filename patterns.

        Returns:
            Dict[str, List[str]]: A dictionary with keys 'R1', 'R2', and 'ignored', each containing lists of filenames.
        """
        if not hasattr(self, "patterns"):
            raise ValueError(
                "No patterns loaded. Either provide a config_path or manually inject patterns."
            )

        # Compile regex patterns from the YAML configuration
        r1_patterns = [
            re.compile(f".*({pattern})(\.|\_|\-).*") for pattern in self.patterns["r1"]
        ]
        r2_patterns = [
            re.compile(f".*({pattern})(\.|\_|\-).*") for pattern in self.patterns["r2"]
        ]
        ignore_patterns = [
            re.compile(f".*({pattern}).*") for pattern in self.patterns["ignore"]
        ]

        # Initialize the result dictionary
        categorized_files = {"R1": [], "R2": [], "ignored": []}

        # Categorize each file
        for filename in self.filenames:
            if any(pattern.search(filename) for pattern in ignore_patterns):
                categorized_files["ignored"].append(filename)
                logger.debug(f"Ignored file: {filename}")
            elif any(pattern.search(filename) for pattern in r1_patterns):
                categorized_files["R1"].append(filename)
                logger.debug(f"Categorized as R1: {filename}")
            elif any(pattern.search(filename) for pattern in r2_patterns):
                categorized_files["R2"].append(filename)
                logger.debug(f"Categorized as R2: {filename}")
            else:
                # If it doesn't match any of the patterns, categorize as ignored
                categorized_files["ignored"].append(filename)
                logger.debug(
                    f"File did not match any R1 or R2 patterns. Index file? {filename}"
                )

        # Sort the filenames alphabetically in each category
        for category in categorized_files:
            categorized_files[category].sort()

        # Check if the number of R1 and R2 files is balanced
        len_r1 = len(categorized_files.get("R1", []))
        len_r2 = len(categorized_files.get("R2", []))

        if len_r1 != len_r2:
            logger.error(f"Unbalanced categories: R1={len_r1}, R2={len_r2}")
            raise ValueError(f"Unbalanced categories: R1={len_r1}, R2={len_r2}")

        return categorized_files
