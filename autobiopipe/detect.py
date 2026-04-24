"""Sequence file detection utilities."""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional, Tuple

FASTA_EXTENSIONS = {".fa", ".fasta", ".fna", ".ffn", ".faa"}
FASTQ_EXTENSIONS = {".fq", ".fastq"}
GZIP_EXTENSIONS = {".gz", ".gzip"}
MATE_PATTERNS: tuple[tuple[re.Pattern[str], str], ...] = (
    (
        re.compile(r"(?i)(^|[._-])R1(?=($|[._-]))"),
        "R2",
    ),
    (
        re.compile(r"(?i)(^|[._-])R2(?=($|[._-]))"),
        "R1",
    ),
    (
        re.compile(r"(?i)(^|[_-])1(?=($|[_-]))"),
        "2",
    ),
    (
        re.compile(r"(?i)(^|[_-])2(?=($|[_-]))"),
        "1",
    ),
)


@dataclass
class DetectionResult:
    """Outcome of file format and layout detection."""

    file_path: Path
    file_type: str
    is_compressed: bool
    sequencing_type: str
    paired_mate: Optional[Path]
    detection_method: str
    confidence: str
    warnings: list[str] = field(default_factory=list)

    def summary(self) -> str:
        """Return a concise, human-readable detection summary.

        Returns:
            A short summary string suitable for logs.
        """
        mate = f", mate={self.paired_mate.name}" if self.paired_mate else ""
        return (
            f"{self.file_path.name}: type={self.file_type}, compressed={self.is_compressed}, "
            f"sequencing={self.sequencing_type}, confidence={self.confidence}{mate}"
        )


def _candidate_extensions(path: Path) -> Tuple[bool, list[str]]:
    """Return compression status and content extensions for a path.

    Args:
        path: File path to inspect.

    Returns:
        A tuple of ``(is_compressed, suffixes_without_gzip)``.
    """
    suffixes = [suffix.lower() for suffix in path.suffixes]
    is_compressed = bool(suffixes and suffixes[-1] in GZIP_EXTENSIONS)
    content_suffixes = suffixes[:-1] if is_compressed else suffixes
    return is_compressed, content_suffixes


def _open_text(path: Path) -> Iterable[str]:
    """Open plain text or gzipped files transparently.

    Args:
        path: File path to read.

    Yields:
        Decoded text lines.
    """
    opener = gzip.open if _candidate_extensions(path)[0] else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        yield from handle


def _first_non_empty_character(path: Path) -> str:
    """Read the first non-whitespace character from a file.

    Args:
        path: File path to inspect.

    Returns:
        The first non-empty character, or an empty string if none is found.
    """
    for line in _open_text(path):
        stripped = line.strip()
        if stripped:
            return stripped[0]
    return ""


def _detect_by_content(path: Path) -> tuple[str, str, str]:
    """Detect format by inspecting the first non-empty character.

    Args:
        path: File path to inspect.

    Returns:
        A tuple of ``(file_type, method, confidence)``.
    """
    first_char = _first_non_empty_character(path)
    if first_char == ">":
        return "fasta", "content-sniff", "high"
    if first_char == "@":
        return "fastq", "content-sniff", "high"
    if not first_char:
        return "unknown", "content-empty", "low"
    return "unknown", "content-ambiguous", "low"


def _detect_by_extension(path: Path) -> tuple[str, str, str]:
    """Detect file type from extensions.

    Args:
        path: File path to inspect.

    Returns:
        A tuple of ``(file_type, method, confidence)``.
    """
    _, suffixes = _candidate_extensions(path)
    ext = suffixes[-1] if suffixes else ""
    if ext in FASTA_EXTENSIONS:
        return "fasta", "extension", "medium"
    if ext in FASTQ_EXTENSIONS:
        return "fastq", "extension", "medium"
    return "unknown", "extension-unknown", "low"


def _candidate_mate_names(name: str) -> list[str]:
    """Generate likely paired-end mate filenames.

    Args:
        name: Input filename.

    Returns:
        Candidate mate filenames inferred from common read pair markers.
    """
    candidates: list[str] = []
    for pattern, replacement in MATE_PATTERNS:
        for match in pattern.finditer(name):
            token_start = match.start() + len(match.group(1))
            token_end = match.end()
            candidates.append(name[:token_start] + replacement + name[token_end:])
    return list(dict.fromkeys(candidates))


def _detect_paired_end(path: Path) -> tuple[str, Optional[Path], list[str]]:
    """Infer paired-end status from filename patterns.

    Args:
        path: File path to inspect.

    Returns:
        A tuple of ``(sequencing_type, mate_path, warnings)``.
    """
    warnings: list[str] = []
    for candidate_name in _candidate_mate_names(path.name):
        candidate_path = path.with_name(candidate_name)
        if candidate_path.exists():
            return "paired-end", candidate_path, warnings
        warnings.append(
            f"Filename suggests paired-end data but mate file '{candidate_name}' was not found."
        )
        return "paired-end", None, warnings
    return "single-end", None, warnings


def detect_file(path: Path) -> DetectionResult:
    """Detect file format, compression, and likely sequencing layout.

    Args:
        path: Path to a candidate FASTA or FASTQ file.

    Returns:
        DetectionResult describing the inferred file properties.

    Raises:
        FileNotFoundError: If the path does not exist.
    """
    resolved_path = Path(path).expanduser().resolve()
    if not resolved_path.exists():
        raise FileNotFoundError(f"File not found: {resolved_path}")

    warnings: list[str] = []
    is_compressed, _ = _candidate_extensions(resolved_path)
    file_type, method, confidence = _detect_by_content(resolved_path)

    ext_type, ext_method, ext_confidence = _detect_by_extension(resolved_path)
    if file_type == "unknown":
        file_type, method, confidence = ext_type, ext_method, ext_confidence
        if file_type == "unknown":
            warnings.append(
                "Unable to identify file type from the first record or file extension."
            )
    elif ext_type not in {"unknown", file_type}:
        warnings.append(
            f"Content suggests '{file_type}' but extension suggests '{ext_type}'. "
            "Using the content-based classification."
        )
        confidence = "medium"

    sequencing_type, paired_mate, pair_warnings = _detect_paired_end(resolved_path)
    warnings.extend(pair_warnings)

    if method == "content-empty":
        warnings.append("File is empty or contains only whitespace.")

    return DetectionResult(
        file_path=resolved_path,
        file_type=file_type,
        is_compressed=is_compressed,
        sequencing_type=sequencing_type,
        paired_mate=paired_mate,
        detection_method=method,
        confidence=confidence,
        warnings=warnings,
    )
