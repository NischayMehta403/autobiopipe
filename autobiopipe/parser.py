"""Memory-efficient FASTA and FASTQ parsing primitives."""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Generator, Iterable, Optional, Union

GZIP_EXTENSIONS = {".gz", ".gzip"}


class ParseError(Exception):
    """Raised when a sequence file is syntactically malformed."""

    def __init__(self, message: str, line_number: Optional[int] = None) -> None:
        """Initialize a ParseError.

        Args:
            message: Human-readable parse error description.
            line_number: Optional 1-based line number where the issue occurred.
        """
        self.line_number = line_number
        location = f" (line {line_number})" if line_number is not None else ""
        super().__init__(f"{message}{location}")


@dataclass
class FastaRecord:
    """A single FASTA record."""

    identifier: str
    description: str
    sequence: str

    def __len__(self) -> int:
        """Return the sequence length."""
        return len(self.sequence)

    def gc_content(self) -> float:
        """Calculate GC content as a percentage.

        Returns:
            GC percentage in the inclusive range 0 to 100.
        """
        if not self.sequence:
            return 0.0
        gc_count = sum(1 for base in self.sequence.upper() if base in {"G", "C"})
        return (gc_count / len(self.sequence)) * 100.0


@dataclass
class FastqRecord:
    """A single FASTQ record."""

    identifier: str
    description: str
    sequence: str
    quality: str

    def __len__(self) -> int:
        """Return the read length."""
        return len(self.sequence)

    @property
    def quality_scores(self) -> list[int]:
        """Decode ASCII-33 quality scores.

        Returns:
            Per-base Phred scores.
        """
        return [ord(symbol) - 33 for symbol in self.quality]

    def avg_quality(self) -> float:
        """Calculate the mean per-base Phred quality score.

        Returns:
            Average decoded quality score, or 0.0 for empty quality strings.
        """
        scores = self.quality_scores
        return sum(scores) / len(scores) if scores else 0.0

    def gc_content(self) -> float:
        """Calculate GC content as a percentage.

        Returns:
            GC percentage in the inclusive range 0 to 100.
        """
        if not self.sequence:
            return 0.0
        gc_count = sum(1 for base in self.sequence.upper() if base in {"G", "C"})
        return (gc_count / len(self.sequence)) * 100.0


AnyRecord = Union[FastaRecord, FastqRecord]


def _open_text(path: Path) -> Iterable[str]:
    """Yield lines from a plain text or gzipped sequence file.

    Args:
        path: File path to read.

    Yields:
        Lines without trailing newline characters.
    """
    opener = gzip.open if path.suffix.lower() in GZIP_EXTENSIONS else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            yield line.rstrip("\r\n")


def _split_header(header: str) -> tuple[str, str]:
    """Split a FASTA or FASTQ header into identifier and description.

    Args:
        header: Raw header text without the leading marker.

    Returns:
        A tuple of ``(identifier, description)``.
    """
    parts = header.strip().split(maxsplit=1)
    identifier = parts[0]
    description = parts[1] if len(parts) > 1 else ""
    return identifier, description


def parse_fasta(
    path: Path,
    max_records: Optional[int] = None,
) -> Generator[FastaRecord, None, None]:
    """Parse FASTA records lazily.

    Args:
        path: FASTA file path.
        max_records: Optional maximum number of records to emit.

    Yields:
        FastaRecord instances.

    Raises:
        FileNotFoundError: If the input file does not exist.
        ParseError: If the FASTA structure is invalid.
    """
    fasta_path = Path(path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    current_header: Optional[str] = None
    sequence_parts: list[str] = []
    header_line: Optional[int] = None
    emitted = 0
    line_number = 0

    for line_number, raw_line in enumerate(_open_text(fasta_path), start=1):
        line = raw_line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if current_header is not None:
                sequence = "".join(sequence_parts).upper()
                if not sequence:
                    raise ParseError(
                        f"Empty sequence for record '{current_header}'",
                        header_line,
                    )
                identifier, description = _split_header(current_header)
                yield FastaRecord(identifier=identifier, description=description, sequence=sequence)
                emitted += 1
                if max_records is not None and emitted >= max_records:
                    return
            current_header = line[1:].strip()
            header_line = line_number
            if not current_header:
                raise ParseError("Empty FASTA header", line_number)
            sequence_parts = []
            continue

        if current_header is None:
            raise ParseError("Sequence encountered before the first FASTA header", line_number)
        sequence_parts.append(line)

    if line_number == 0:
        raise ParseError("FASTA file is empty")
    if current_header is None:
        raise ParseError("No FASTA header found")

    sequence = "".join(sequence_parts).upper()
    if not sequence:
        raise ParseError(f"Empty sequence for record '{current_header}'", header_line)
    identifier, description = _split_header(current_header)
    yield FastaRecord(identifier=identifier, description=description, sequence=sequence)


def parse_fastq(
    path: Path,
    max_records: Optional[int] = None,
) -> Generator[FastqRecord, None, None]:
    """Parse FASTQ records lazily.

    Args:
        path: FASTQ file path.
        max_records: Optional maximum number of records to emit.

    Yields:
        FastqRecord instances.

    Raises:
        FileNotFoundError: If the input file does not exist.
        ParseError: If the FASTQ structure is invalid.
    """
    fastq_path = Path(path)
    if not fastq_path.exists():
        raise FileNotFoundError(f"FASTQ file not found: {fastq_path}")

    lines = iter(enumerate(_open_text(fastq_path), start=1))
    emitted = 0
    saw_content = False

    while True:
        try:
            header_line_number, header_line = next(lines)
        except StopIteration:
            if not saw_content:
                raise ParseError("FASTQ file is empty")
            return

        if not header_line.strip():
            continue

        saw_content = True
        if not header_line.startswith("@"):
            raise ParseError("FASTQ record header must start with '@'", header_line_number)

        header_text = header_line[1:].strip()
        if not header_text:
            raise ParseError("FASTQ header is empty", header_line_number)

        try:
            sequence_line_number, sequence_line = next(lines)
            separator_line_number, separator_line = next(lines)
            quality_line_number, quality_line = next(lines)
        except StopIteration as exc:
            raise ParseError(
                "FASTQ record is incomplete; expected four lines per record",
                header_line_number,
            ) from exc

        sequence = sequence_line.strip().upper()
        if not sequence:
            raise ParseError("FASTQ sequence line is empty", sequence_line_number)
        if not separator_line.startswith("+"):
            raise ParseError(
                "FASTQ separator line must start with '+'",
                separator_line_number,
            )

        quality = quality_line.strip()
        if len(sequence) != len(quality):
            raise ParseError(
                "FASTQ quality length does not match sequence length",
                quality_line_number,
            )

        identifier, description = _split_header(header_text)
        yield FastqRecord(
            identifier=identifier,
            description=description,
            sequence=sequence,
            quality=quality,
        )
        emitted += 1
        if max_records is not None and emitted >= max_records:
            return


def parse_file(
    path: Path,
    file_type: str,
    max_records: Optional[int] = None,
) -> Generator[AnyRecord, None, None]:
    """Dispatch to the appropriate parser based on file type.

    Args:
        path: Input file path.
        file_type: File type label, typically ``fasta`` or ``fastq``.
        max_records: Optional cap on emitted records.

    Yields:
        Parsed FastaRecord or FastqRecord instances.

    Raises:
        ValueError: If the requested file type is unsupported.
    """
    normalized = file_type.lower()
    if normalized == "fasta":
        yield from parse_fasta(path, max_records=max_records)
        return
    if normalized == "fastq":
        yield from parse_fastq(path, max_records=max_records)
        return
    raise ValueError(f"Unsupported file type: {file_type}")
