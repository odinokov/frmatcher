"""Regression tests locking frmatcher's public behavior.

These pin the documented CLI/API contract so future edits can be checked
against current behavior.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from frmatcher import (
    FASTQ_EXTENSIONS,
    FastqPairingError,
    PairingStats,
    group_fastq_pairs,
    is_index_file,
    is_index_read_name,
    iter_fastq_files,
    parse_fastq_name,
    write_fastq_pairs,
    __version__,
)


def _touch(directory: Path, *names: str) -> None:
    for name in names:
        (directory / name).touch()


# --------------------------------------------------------------------------- #
# parse_fastq_name
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "filename, read, sep",
    [
        ("sample_R1_001.fastq.gz", "1", "_"),
        ("sample_R2_001.fastq.gz", "2", "_"),
        ("sample_1.fastq.gz", "1", "_"),
        ("sample_2.fq", "2", "_"),
        ("sample.R1.fastq", "1", "."),
        ("sample-R2.fq.bz2", "2", "-"),
    ],
)
def test_parse_fastq_name_reads(filename, read, sep):
    parsed = parse_fastq_name(filename)
    assert parsed is not None
    assert parsed["read"] == read
    assert parsed["sep"] == sep


def test_parse_fastq_name_greedy_prefix_last_read_wins():
    # Documented behavior: greedy prefix => last _R1/_1 wins.
    parsed = parse_fastq_name("sample_1_R2.fastq.gz")
    assert parsed is not None
    assert parsed["read"] == "2"
    assert parsed["prefix"] == "sample_1"


@pytest.mark.parametrize(
    "filename",
    ["notes.txt", "sample.fastq.gz", "sample_R3_001.fastq.gz", "README.md"],
)
def test_parse_fastq_name_non_matches(filename):
    assert parse_fastq_name(filename) is None


# --------------------------------------------------------------------------- #
# index detection
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "filename, expected",
    [
        ("sample_I1_001.fastq.gz", True),
        ("sample_I2.fastq.gz", True),
        ("sample_R1_001.fastq.gz", False),
        ("sample_1.fastq.gz", False),
        ("patient_I1_R1.fastq.gz", False),
    ],
)
def test_is_index_read_name(filename, expected):
    assert is_index_read_name(filename) is expected


@pytest.mark.parametrize(
    "stemish, expected",
    [("sample_I1", True), ("I2_sample", True), ("sample_R1", False), ("plain", False)],
)
def test_is_index_file(stemish, expected):
    assert is_index_file(stemish) is expected


# --------------------------------------------------------------------------- #
# iter_fastq_files
# --------------------------------------------------------------------------- #
def test_iter_fastq_files_extensions_and_case(tmp_path):
    _touch(
        tmp_path,
        "a_R1.fastq.gz",
        "b_R1.FASTQ.GZ",  # case-insensitive
        "c_R1.fq",
        "ignore.txt",
        "ignore.bam",
    )
    found = {p.name for p in iter_fastq_files(tmp_path)}
    assert found == {"a_R1.fastq.gz", "b_R1.FASTQ.GZ", "c_R1.fq"}


def test_iter_fastq_files_non_recursive_skips_subdirs(tmp_path):
    _touch(tmp_path, "top_R1.fastq.gz")
    sub = tmp_path / "sub"
    sub.mkdir()
    _touch(sub, "deep_R1.fastq.gz")

    flat = {p.name for p in iter_fastq_files(tmp_path)}
    deep = {p.name for p in iter_fastq_files(tmp_path, recursive=True)}
    assert flat == {"top_R1.fastq.gz"}
    assert deep == {"top_R1.fastq.gz", "deep_R1.fastq.gz"}


def test_iter_fastq_files_custom_extensions(tmp_path):
    _touch(tmp_path, "x_R1.fq.zst", "y_R1.fastq.gz")
    custom = FASTQ_EXTENSIONS | frozenset({".fq.zst"})
    found = {p.name for p in iter_fastq_files(tmp_path, extensions=custom)}
    assert found == {"x_R1.fq.zst", "y_R1.fastq.gz"}


# --------------------------------------------------------------------------- #
# group_fastq_pairs
# --------------------------------------------------------------------------- #
def test_group_fastq_pairs_counts_and_grouping(tmp_path):
    _touch(
        tmp_path,
        "sampleA_R1_001.fastq.gz",
        "sampleA_R2_001.fastq.gz",
        "sampleB_I1_001.fastq.gz",  # index -> skipped
        "sampleB_R1_001.fastq.gz",
        "sampleB_R2_001.fastq.gz",
        "weird.fastq.gz",  # unmatched
    )
    stats = group_fastq_pairs(tmp_path)
    assert isinstance(stats, PairingStats)
    assert stats.n_fastq == 6
    assert stats.n_index == 1
    assert stats.n_unmatched == 1
    # Two paired groups; R1 and R2 of a sample share one key.
    assert len(stats.groups) == 2
    for reads in stats.groups.values():
        assert len(reads["1"]) == 1
        assert len(reads["2"]) == 1


def test_group_fastq_pairs_recursive_keeps_directories_separate(tmp_path):
    a = tmp_path / "a"
    b = tmp_path / "b"
    a.mkdir()
    b.mkdir()
    _touch(a, "sample_R1.fastq.gz")
    _touch(b, "sample_R2.fastq.gz")

    stats = group_fastq_pairs(tmp_path, recursive=True)

    assert len(stats.groups) == 2
    assert all(len(reads["1"]) != len(reads["2"]) for reads in stats.groups.values())


def test_group_fastq_pairs_sample_i1_token_is_not_index_read(tmp_path):
    _touch(tmp_path, "patient_I1_R1.fastq.gz", "patient_I1_R2.fastq.gz")

    stats = group_fastq_pairs(tmp_path)

    assert stats.n_index == 0
    assert len(stats.groups) == 1
    reads = next(iter(stats.groups.values()))
    assert len(reads["1"]) == 1
    assert len(reads["2"]) == 1


# --------------------------------------------------------------------------- #
# write_fastq_pairs
# --------------------------------------------------------------------------- #
def test_write_fastq_pairs_basic(tmp_path):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    _touch(in_dir, "s_R1_001.fastq.gz", "s_R2_001.fastq.gz")
    out = tmp_path / "out" / "pairs.tsv"  # parent does not exist yet

    write_fastq_pairs(in_dir, out)

    assert out.exists()
    assert not out.with_name(out.name + ".tmp").exists()  # atomic: no temp left
    rows = out.read_text().splitlines()
    assert len(rows) == 1
    r1, r2 = rows[0].split("\t")
    assert r1.endswith("s_R1_001.fastq.gz")
    assert r2.endswith("s_R2_001.fastq.gz")


def test_write_fastq_pairs_not_a_directory(tmp_path):
    with pytest.raises(NotADirectoryError):
        write_fastq_pairs(tmp_path / "does_not_exist", tmp_path / "out.tsv")


def test_write_fastq_pairs_strict_unpaired(tmp_path):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    _touch(in_dir, "s_R1_001.fastq.gz")  # R1 only -> unpaired
    out = tmp_path / "pairs.tsv"

    with pytest.raises(FastqPairingError):
        write_fastq_pairs(in_dir, out, strict=True)
    assert not out.with_name(out.name + ".tmp").exists()


def test_write_fastq_pairs_strict_unmatched(tmp_path):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    _touch(in_dir, "p_R1.fastq.gz", "p_R2.fastq.gz", "weird.fastq.gz")
    out = tmp_path / "pairs.tsv"

    with pytest.raises(FastqPairingError):
        write_fastq_pairs(in_dir, out, strict=True)


def test_write_fastq_pairs_relative_paths_preserved(tmp_path, monkeypatch):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    _touch(in_dir, "s_R1.fastq.gz", "s_R2.fastq.gz")
    monkeypatch.chdir(tmp_path)

    write_fastq_pairs("in", "pairs.tsv")  # relative in-dir
    rows = (tmp_path / "pairs.tsv").read_text().splitlines()
    assert rows[0].split("\t")[0].startswith("in/")  # relative preserved


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def test_cli_version():
    out = subprocess.run(
        [sys.executable, "-m", "frmatcher", "--version"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert out.stdout.strip() == f"frmatcher {__version__}"


def test_cli_end_to_end(tmp_path):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    _touch(in_dir, "s_R1_001.fastq.gz", "s_R2_001.fastq.gz")
    out = tmp_path / "pairs.tsv"

    subprocess.run(
        [sys.executable, "-m", "frmatcher", "-i", str(in_dir), "-o", str(out)],
        capture_output=True,
        text=True,
        check=True,
    )
    rows = out.read_text().splitlines()
    assert len(rows) == 1
    assert rows[0].split("\t")[0].endswith("s_R1_001.fastq.gz")


def test_cli_strict_error_without_traceback(tmp_path):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    _touch(in_dir, "s_R1.fastq.gz")
    out = tmp_path / "pairs.tsv"

    result = subprocess.run(
        [sys.executable, "-m", "frmatcher", "-i", str(in_dir), "-o", str(out), "--strict"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 1
    assert "unpaired group" in result.stderr
    assert "Traceback" not in result.stderr
