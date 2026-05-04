"""Pair FASTQ R1/R2 files and output as TSV."""
from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterator, NamedTuple, Optional, Union

from loguru import logger

logger.disable(__name__)

__version__ = "0.1.0"

FASTQ_EXTENSIONS: frozenset[str] = frozenset(
    {".fq", ".fastq", ".fq.gz", ".fastq.gz", ".fq.bz2", ".fastq.bz2"}
)

_READ_PATTERN = re.compile(
    r"^(?P<prefix>.+)(?P<sep>[._-])R?(?P<read>[12])(?P<tail>(?:[._-].*)?)"
    r"\.(?:f(?:ast)?q)(?:\.(?:gz|bz2))?$",
    re.IGNORECASE,
)

_INDEX_READ_PATTERN = re.compile(
    r"^.+[._-]I[12](?:[._-].*)?\.(?:f(?:ast)?q)(?:\.(?:gz|bz2))?$",
    re.IGNORECASE,
)

_INDEX_PATTERN = re.compile(r"(?:^|[._-])I[12](?:[._-]|$)", re.IGNORECASE)


class FastqPairingError(ValueError):
    pass


class PairingStats(NamedTuple):
    n_fastq: int
    n_index: int
    n_unmatched: int
    groups: dict[str, dict[str, list[str]]]


def iter_fastq_files(
    directory: Path,
    recursive: bool = False,
    extensions: frozenset[str] = FASTQ_EXTENSIONS,
) -> Iterator[Path]:
    iterator = directory.rglob("*") if recursive else directory.iterdir()
    for path in iterator:
        if path.is_file():
            name = path.name.lower()
            if any(name.endswith(ext) for ext in extensions):
                yield path


def parse_fastq_name(filename: str) -> Optional[dict[str, str]]:
    # greedy prefix => last _R1/_1 wins (e.g., sample_1_R2.fastq.gz parses as R2)
    m = _READ_PATTERN.match(filename)
    return m.groupdict() if m else None


def is_index_file(stemish: str) -> bool:
    return bool(_INDEX_PATTERN.search(stemish))


def is_index_read_name(filename: str) -> bool:
    return bool(_INDEX_READ_PATTERN.match(filename))


def group_fastq_pairs(directory: Path, recursive: bool = False) -> PairingStats:
    groups: dict[str, dict[str, list[str]]] = defaultdict(lambda: {"1": [], "2": []})
    n_fastq = n_index = n_unmatched = 0

    for path in iter_fastq_files(directory, recursive):
        n_fastq += 1
        if is_index_read_name(path.name):
            n_index += 1
            continue

        parsed = parse_fastq_name(path.name)
        if not parsed:
            n_unmatched += 1
            logger.debug(f"Unmatched: {path.name}")
            continue
        stemish = f"{parsed['prefix']}{parsed['tail']}"
        if is_index_file(stemish):
            n_index += 1
            continue
        key = f"{parsed['prefix']}{parsed['sep']}{{read}}{parsed['tail']}"
        groups[key][parsed["read"]].append(str(path))

    return PairingStats(n_fastq=n_fastq, n_index=n_index, n_unmatched=n_unmatched, groups=dict(groups))


def _pair_rows(
    groups: dict[str, dict[str, list[str]]],
) -> tuple[list[tuple[str, str]], list[tuple[str, int, int]]]:
    pairs: list[tuple[str, str]] = []
    unpaired: list[tuple[str, int, int]] = []
    for key in sorted(groups):
        r1 = sorted(groups[key]["1"])
        r2 = sorted(groups[key]["2"])
        if not r1 or len(r1) != len(r2):
            unpaired.append((key, len(r1), len(r2)))
        else:
            pairs.extend(zip(r1, r2))
    return pairs, unpaired


def write_fastq_pairs(
    in_dir: Union[str, Path],
    out_tsv: Union[str, Path],
    recursive: bool = False,
    strict: bool = False,
) -> None:
    in_path = Path(in_dir)
    out_path = Path(out_tsv)

    if not in_path.is_dir():
        raise NotADirectoryError(f"Not a directory: {in_dir}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"Scanning {in_path} (recursive={recursive})")

    stats = group_fastq_pairs(in_path, recursive)
    pairs, unpaired = _pair_rows(stats.groups)

    tmp = out_path.with_name(out_path.name + ".tmp")
    try:
        with tmp.open("w", encoding="utf-8") as fh:
            for a, b in pairs:
                fh.write(f"{a}\t{b}\n")
    except Exception:
        tmp.unlink(missing_ok=True)
        raise

    if unpaired:
        shown = unpaired[:20]
        extra = f" +{len(unpaired) - 20} more" if len(unpaired) > 20 else ""
        logger.warning(f"Unpaired groups (key, nR1, nR2): {shown}{extra}")
        if strict:
            tmp.unlink(missing_ok=True)
            raise FastqPairingError(f"{len(unpaired)} unpaired group(s); first: {unpaired[0]}")

    if stats.n_unmatched and strict:
        tmp.unlink(missing_ok=True)
        raise FastqPairingError(f"{stats.n_unmatched} unmatched FASTQ file(s)")

    tmp.replace(out_path)
    logger.info(
        f"Done: fastq={stats.n_fastq}, unmatched={stats.n_unmatched}, "
        f"index={stats.n_index}, pairs={len(pairs)}, unpaired={len(unpaired)}"
    )


def main() -> int:
    p = argparse.ArgumentParser(prog="frmatcher", description="Pair FASTQ R1/R2 files and write TSV (R1<TAB>R2).")
    p.add_argument("-i", "--in-dir", required=True, help="Directory containing FASTQ files")
    p.add_argument("-o", "--out-tsv", required=True, help="Output TSV file path")
    p.add_argument("-r", "--recursive", action="store_true", help="Search subdirectories")
    p.add_argument("--strict", action="store_true", help="Fail if any unpaired group exists")
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=["TRACE", "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    )
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = p.parse_args()

    logger.enable(__name__)
    logger.remove()
    logger.add(sys.stderr, level=args.log_level)

    write_fastq_pairs(args.in_dir, args.out_tsv, recursive=args.recursive, strict=args.strict)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
