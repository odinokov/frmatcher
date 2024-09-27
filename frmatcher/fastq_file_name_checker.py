import os
import re
from typing import Dict, List, Optional

import yaml
from loguru import logger


class FilenameLengthMismatchError(ValueError):
    """Custom exception for filename length mismatches."""


class PatternsNotLoadedError(ValueError):
    """Custom exception when patterns are not loaded."""


class FastqFileNameChecker:
    def __init__(
        self,
        filenames: List[str],
        config_path: Optional[str] = None,
        length_check: bool = False,
        verbose: bool = False,
    ) -> None:
        """
        Initialize the FastqFileNameChecker with a list of filenames.

        Args:
            filenames (List[str]): List of file paths to categorize.
            config_path (str, optional): Path to the YAML configuration file.
                Defaults to None, and patterns must be manually injected.
            length_check (bool): Whether to check if all filenames have the same length.
                Default is False.
            verbose (bool): Whether to enable detailed logging.
                Default is False.

        Raises:
            FilenameLengthMismatchError: If filenames have different lengths and length_check is True.
            FileNotFoundError: If the configuration file does not exist.
            yaml.YAMLError: If the configuration file is invalid.
        """
        self.filenames = filenames
        self.verbose = verbose

        # Configure logging
        logger.remove()  # Remove the default handler
        if self.verbose:
            logger.add(lambda msg: print(msg, end=""), level="DEBUG", colorize=True)
        else:
            logger.add(lambda msg: print(msg, end=""), level="ERROR", colorize=True)

        # Only load patterns from config if no patterns are directly injected
        if config_path is None and not hasattr(self, "patterns"):
            self.patterns = None  # Set patterns to None initially
        elif config_path is not None:
            self.patterns = self.load_patterns(config_path)

        # Precompile regex patterns if patterns are loaded or injected
        if self.patterns:
            self.compiled_patterns = self.compile_patterns(self.patterns)

        if length_check:
            self._check_filename_lengths()

    def load_patterns(self, config_path: str) -> Dict[str, List[str]]:
        """
        Load patterns from a YAML configuration file.

        Args:
            config_path (str): Path to the YAML configuration file.

        Returns:
            Dict[str, List[str]]: A dictionary with R1, R2, and ignore patterns.

        Raises:
            FileNotFoundError: If the configuration file does not exist.
            yaml.YAMLError: If the configuration file is invalid.
        """
        if not os.path.isfile(config_path):
            logger.error(f"Configuration file not found at {config_path}")
            raise FileNotFoundError(f"Configuration file not found at {config_path}")

        try:
            with open(config_path, "r") as file:
                config = yaml.safe_load(file)
            return config["patterns"]
        except yaml.YAMLError as e:
            logger.error(f"Error parsing YAML configuration: {e}")
            raise

    def compile_patterns(
        self, patterns: Dict[str, List[str]]
    ) -> Dict[str, List[re.Pattern]]:
        """
        Compile regex patterns for R1, R2, and ignore categories.

        Args:
            patterns (Dict[str, List[str]]): Raw patterns from configuration.

        Returns:
            Dict[str, List[re.Pattern]]: Compiled regex patterns.
        """
        compiled = {
            "R1": [
                re.compile(f".*({pattern})([._-]).*")
                for pattern in patterns.get("r1", [])
            ],
            "R2": [
                re.compile(f".*({pattern})([._-]).*")
                for pattern in patterns.get("r2", [])
            ],
            "ignore": [
                re.compile(f".*({pattern}).*") for pattern in patterns.get("ignore", [])
            ],
        }
        return compiled

    def _check_filename_lengths(self) -> None:
        """
        Checks if all filenames have the same length.

        Raises:
            FilenameLengthMismatchError: If filenames do not have the same length.
        """
        lengths = list(map(len, self.filenames))
        if len(set(lengths)) > 1:
            logger.error(
                "Filenames do not all have the same length. Please ensure all filenames are consistent."
            )
            raise FilenameLengthMismatchError(
                "Filenames do not all have the same length. Please ensure all filenames are consistent."
            )

    def categorize_fastq_files(self) -> Dict[str, List[str]]:
        """
        Categorizes FASTQ files into R1, R2, or ignored based on filename patterns.

        Returns:
            Dict[str, List[str]]: A dictionary with keys 'R1', 'R2', and 'ignored',
                each containing lists of full file paths.

        Raises:
            PatternsNotLoadedError: If regex patterns are not loaded.
            FilenameLengthMismatchError: If the number of R1 and R2 files is unbalanced.
        """
        if not hasattr(self, "compiled_patterns"):
            raise PatternsNotLoadedError("Patterns not loaded or compiled.")

        categorized_files: Dict[str, List[str]] = {"R1": [], "R2": [], "ignored": []}

        for full_path in self.filenames:
            filename = os.path.basename(full_path)
            if any(
                pattern.search(filename) for pattern in self.compiled_patterns["ignore"]
            ):
                categorized_files["ignored"].append(full_path)
            elif any(
                pattern.search(filename) for pattern in self.compiled_patterns["R1"]
            ):
                categorized_files["R1"].append(full_path)
            elif any(
                pattern.search(filename) for pattern in self.compiled_patterns["R2"]
            ):
                categorized_files["R2"].append(full_path)
            else:
                categorized_files["ignored"].append(full_path)

        # Sort the filenames alphabetically in each category
        for category in categorized_files:
            categorized_files[category].sort()

        # Check if the number of R1 and R2 files is balanced
        len_r1 = len(categorized_files.get("R1", []))
        len_r2 = len(categorized_files.get("R2", []))

        if len_r1 != len_r2:
            logger.error(f"Unbalanced categories: R1={len_r1}, R2={len_r2}")
            raise FilenameLengthMismatchError(
                f"Unbalanced categories: R1={len_r1}, R2={len_r2}"
            )

        return categorized_files
