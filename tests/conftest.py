"""Shared fixtures for all AutoBioPipe tests."""

import textwrap
from pathlib import Path
import pytest


@pytest.fixture()
def fasta_file(tmp_path: Path) -> Path:
    """A small valid FASTA file."""
    content = textwrap.dedent("""\
        >seq1 human chromosome 1
        ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
        >seq2 gc balanced
        GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
        >seq3 at rich
        ATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT
    """)
    p = tmp_path / "sample.fasta"
    p.write_text(content)
    return p

@pytest.fixture()
def fastq_file(tmp_path: Path) -> Path:
    """A small valid FASTQ file."""
    p = tmp_path / "sample.fastq"
    p.write_text(
        "@read1\n"
        "ATGC\n"
        "+\n"
        "IIII\n"
        "@read2\n"
        "GCGC\n"
        "+\n"
        "5555\n"
        "@read3\n"
        "ATAT\n"
        "+\n"
        "!!!!\n"
    )
    return p

@pytest.fixture()
def malformed_fastq_file(tmp_path: Path) -> Path:
    """FASTQ with mismatched quality length."""
    content = textwrap.dedent("""\
        @bad_read
        ATGCATGC
        +
        III
    """)
    p = tmp_path / "malformed.fastq"
    p.write_text(content)
    return p


@pytest.fixture()
def empty_file(tmp_path: Path) -> Path:
    """An empty file."""
    p = tmp_path / "empty.fasta"
    p.write_text("")
    return p
