"""Tests for autobiopipe.parser"""

from pathlib import Path
import pytest

from autobiopipe.parser import (
    FastaRecord,
    FastqRecord,
    ParseError,
    parse_fasta,
    parse_fastq,
    parse_file,
)


# ── FastaRecord ──────────────────────────────────────────────


class TestFastaRecord:
    def test_len(self) -> None:
        rec = FastaRecord("s1", "s1", "ATGCATGC")
        assert len(rec) == 8

    def test_gc_pure_gc(self) -> None:
        rec = FastaRecord("s", "s", "GCGCGC")
        assert rec.gc_content() == pytest.approx(100.0)

    def test_gc_pure_at(self) -> None:
        rec = FastaRecord("s", "s", "ATATAT")
        assert rec.gc_content() == pytest.approx(0.0)

    def test_gc_balanced(self) -> None:
        rec = FastaRecord("s", "s", "ATGC")
        assert rec.gc_content() == pytest.approx(50.0)

    def test_gc_empty(self) -> None:
        rec = FastaRecord("s", "s", "")
        assert rec.gc_content() == pytest.approx(0.0)


# ── FastqRecord ──────────────────────────────────────────────


class TestFastqRecord:
    def test_quality_decoded(self) -> None:
        # 'I' = ASCII 73 → Phred 40
        rec = FastqRecord("r", "r", "ATGC", "IIII")
        assert rec.quality_scores == [40, 40, 40, 40]

    def test_avg_quality_high(self) -> None:
        rec = FastqRecord("r", "r", "ATGC", "IIII")
        assert rec.avg_quality() == pytest.approx(40.0)

    def test_avg_quality_low(self) -> None:
        # '!' = Phred 0
        rec = FastqRecord("r", "r", "ATGC", "!!!!")
        assert rec.avg_quality() == pytest.approx(0.0)

    def test_avg_quality_mixed(self) -> None:
        # '!' = 0, 'I' = 40 → mean = 20
        rec = FastqRecord("r", "r", "AT", "!I")
        assert rec.avg_quality() == pytest.approx(20.0)


# ── parse_fasta ──────────────────────────────────────────────


class TestParseFasta:
    def test_parses_three_records(self, fasta_file: Path) -> None:
        records = list(parse_fasta(fasta_file))
        assert len(records) == 3

    def test_correct_identifiers(self, fasta_file: Path) -> None:
        records = list(parse_fasta(fasta_file))
        assert records[0].identifier == "seq1"
        assert records[1].identifier == "seq2"

    def test_sequence_uppercase(self, fasta_file: Path) -> None:
        records = list(parse_fasta(fasta_file))
        for rec in records:
            assert rec.sequence == rec.sequence.upper()

    def test_max_records_cap(self, fasta_file: Path) -> None:
        records = list(parse_fasta(fasta_file, max_records=2))
        assert len(records) == 2

    def test_file_not_found(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            list(parse_fasta(tmp_path / "missing.fasta"))

    def test_empty_file_raises(self, empty_file: Path) -> None:
        with pytest.raises(ParseError):
            list(parse_fasta(empty_file))

    def test_multiline_sequence_joined(self, tmp_path: Path) -> None:
        p = tmp_path / "multi.fasta"
        p.write_text(">seq1\nATGC\nATGC\n")
        records = list(parse_fasta(p))
        assert records[0].sequence == "ATGCATGC"

    def test_blank_lines_skipped(self, tmp_path: Path) -> None:
        p = tmp_path / "blanks.fasta"
        p.write_text(">seq1\n\nATGC\n\n>seq2\nGCGC\n")
        records = list(parse_fasta(p))
        assert len(records) == 2


# ── parse_fastq ──────────────────────────────────────────────


class TestParseFastq:
    def test_parses_three_records(self, fastq_file: Path) -> None:
        records = list(parse_fastq(fastq_file))
        assert len(records) == 3

    def test_correct_identifier(self, fastq_file: Path) -> None:
        records = list(parse_fastq(fastq_file))
        assert records[0].identifier == "read1"

    def test_quality_decoded(self, fastq_file: Path) -> None:
        records = list(parse_fastq(fastq_file))
        # First read uses 'I' = Phred 40
        assert all(s == 40 for s in records[0].quality_scores)

    def test_max_records_cap(self, fastq_file: Path) -> None:
        records = list(parse_fastq(fastq_file, max_records=1))
        assert len(records) == 1

    def test_file_not_found(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            list(parse_fastq(tmp_path / "missing.fastq"))

    def test_malformed_quality_raises(self, malformed_fastq_file: Path) -> None:
        with pytest.raises(ParseError, match="Quality length"):
            list(parse_fastq(malformed_fastq_file))

    def test_missing_at_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "bad.fastq"
        p.write_text("read1\nATGC\n+\nIIII\n")
        with pytest.raises(ParseError, match="Expected '@'"):
            list(parse_fastq(p))


# ── parse_file dispatcher ─────────────────────────────────────


class TestParseFile:
    def test_dispatches_fasta(self, fasta_file: Path) -> None:
        records = list(parse_file(fasta_file, "fasta"))
        assert all(isinstance(r, FastaRecord) for r in records)

    def test_dispatches_fastq(self, fastq_file: Path) -> None:
        records = list(parse_file(fastq_file, "fastq"))
        assert all(isinstance(r, FastqRecord) for r in records)

    def test_invalid_type_raises(self, fasta_file: Path) -> None:
        with pytest.raises(ValueError, match="Unsupported"):
            list(parse_file(fasta_file, "bam"))
