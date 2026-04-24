"""Shared pytest fixtures for AutoBioPipe."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest


@pytest.fixture()
def fasta_file(tmp_path: Path) -> Path:
    """Create a valid FASTA file with three records."""
    path = tmp_path / "sample.fasta"
    path.write_text(
        ">seq1 alpha\n"
        "ATGC\n"
        ">seq2 beta\n"
        "GGCC\n"
        ">seq3 gamma\n"
        "AATT\n",
        encoding="utf-8",
    )
    return path


@pytest.fixture()
def multiline_fasta_file(tmp_path: Path) -> Path:
    """Create a FASTA file with multiline sequences."""
    path = tmp_path / "multiline.fasta"
    path.write_text(
        ">seq1 long description\n"
        "AT\n"
        "GC\n"
        ">seq2 another description\n"
        "GG\n"
        "CC\n",
        encoding="utf-8",
    )
    return path


@pytest.fixture()
def fastq_file(tmp_path: Path) -> Path:
    """Create a valid FASTQ file with three short reads."""
    path = tmp_path / "sample.fastq"
    path.write_text(
        "@read1 alpha\n"
        "ATGC\n"
        "+\n"
        "IIII\n"
        "@read2 beta\n"
        "GGCC\n"
        "+\n"
        "5555\n"
        "@read3 gamma\n"
        "AATT\n"
        "+\n"
        "!!!!\n",
        encoding="utf-8",
    )
    return path


@pytest.fixture()
def gzip_fastq_file(tmp_path: Path) -> Path:
    """Create a gzipped FASTQ file."""
    path = tmp_path / "sample.fastq.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(
            "@read1 alpha\n"
            "ATGC\n"
            "+\n"
            "IIII\n"
            "@read2 beta\n"
            "GGCC\n"
            "+\n"
            "5555\n"
        )
    return path


@pytest.fixture()
def gzip_fasta_file(tmp_path: Path) -> Path:
    """Create a gzipped FASTA file."""
    path = tmp_path / "sample.fasta.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(">seq1\nATGC\n>seq2\nGGCC\n")
    return path


@pytest.fixture()
def malformed_fastq_file(tmp_path: Path) -> Path:
    """Create a FASTQ file with mismatched quality length."""
    path = tmp_path / "malformed.fastq"
    path.write_text("@bad\nATGC\n+\nIII\n", encoding="utf-8")
    return path


@pytest.fixture()
def incomplete_fastq_file(tmp_path: Path) -> Path:
    """Create a truncated FASTQ file."""
    path = tmp_path / "incomplete.fastq"
    path.write_text("@bad\nATGC\n+\n", encoding="utf-8")
    return path


@pytest.fixture()
def empty_file(tmp_path: Path) -> Path:
    """Create an empty file."""
    path = tmp_path / "empty.fasta"
    path.write_text("", encoding="utf-8")
    return path


@pytest.fixture()
def weird_extension_fastq(tmp_path: Path) -> Path:
    """Create a FASTQ file with an uncommon extension."""
    path = tmp_path / "reads.data"
    path.write_text("@read1\nATGC\n+\nIIII\n", encoding="utf-8")
    return path


@pytest.fixture()
def paired_fastq_files(tmp_path: Path) -> tuple[Path, Path]:
    """Create an R1/R2 paired FASTQ file set."""
    r1 = tmp_path / "sample_R1.fastq"
    r2 = tmp_path / "sample_R2.fastq"
    content = "@read1\nATGC\n+\nIIII\n"
    r1.write_text(content, encoding="utf-8")
    r2.write_text(content, encoding="utf-8")
    return r1, r2
