"""Tests for parsing FASTA and FASTQ inputs."""

from __future__ import annotations

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


class TestFastaRecord:
    def test_length(self) -> None:
        assert len(FastaRecord("id", "desc", "ATGC")) == 4

    def test_gc_content_balanced(self) -> None:
        assert FastaRecord("id", "desc", "ATGC").gc_content() == pytest.approx(50.0)

    def test_gc_content_zero_for_empty(self) -> None:
        assert FastaRecord("id", "desc", "").gc_content() == pytest.approx(0.0)

    def test_gc_content_all_gc(self) -> None:
        assert FastaRecord("id", "desc", "GGCC").gc_content() == pytest.approx(100.0)

    def test_gc_content_all_at(self) -> None:
        assert FastaRecord("id", "desc", "AATT").gc_content() == pytest.approx(0.0)


class TestFastqRecord:
    def test_length(self) -> None:
        assert len(FastqRecord("id", "desc", "ATGC", "IIII")) == 4

    def test_quality_scores_ascii33(self) -> None:
        assert FastqRecord("id", "desc", "ATGC", "IIII").quality_scores == [40, 40, 40, 40]

    def test_average_quality(self) -> None:
        assert FastqRecord("id", "desc", "ATGC", "!I!I").avg_quality() == pytest.approx(20.0)

    def test_gc_content(self) -> None:
        assert FastqRecord("id", "desc", "GGCC", "IIII").gc_content() == pytest.approx(100.0)

    def test_avg_quality_zero_for_empty_quality(self) -> None:
        assert FastqRecord("id", "desc", "", "").avg_quality() == pytest.approx(0.0)


class TestParseFasta:
    def test_parses_three_records(self, fasta_file: Path) -> None:
        assert len(list(parse_fasta(fasta_file))) == 3

    def test_identifier_is_first_token(self, fasta_file: Path) -> None:
        assert list(parse_fasta(fasta_file))[0].identifier == "seq1"

    def test_description_is_rest_of_header(self, fasta_file: Path) -> None:
        assert list(parse_fasta(fasta_file))[0].description == "alpha"

    def test_multiline_sequences_are_joined(self, multiline_fasta_file: Path) -> None:
        assert list(parse_fasta(multiline_fasta_file))[0].sequence == "ATGC"

    def test_sequences_are_uppercased(self, tmp_path: Path) -> None:
        path = tmp_path / "lower.fasta"
        path.write_text(">seq1\natgc\n", encoding="utf-8")
        assert list(parse_fasta(path))[0].sequence == "ATGC"

    def test_max_records_caps_output(self, fasta_file: Path) -> None:
        assert len(list(parse_fasta(fasta_file, max_records=2))) == 2

    def test_blank_lines_are_skipped(self, tmp_path: Path) -> None:
        path = tmp_path / "blank.fasta"
        path.write_text(">seq1\n\nATGC\n\n>seq2\nGGCC\n", encoding="utf-8")
        assert len(list(parse_fasta(path))) == 2

    def test_raises_for_missing_file(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            list(parse_fasta(tmp_path / "missing.fasta"))

    def test_raises_for_empty_file(self, empty_file: Path) -> None:
        with pytest.raises(ParseError, match="empty"):
            list(parse_fasta(empty_file))

    def test_raises_for_missing_header(self, tmp_path: Path) -> None:
        path = tmp_path / "bad.fasta"
        path.write_text("ATGC\n", encoding="utf-8")
        with pytest.raises(ParseError, match="first FASTA header"):
            list(parse_fasta(path))

    def test_raises_for_empty_header(self, tmp_path: Path) -> None:
        path = tmp_path / "bad.fasta"
        path.write_text(">\nATGC\n", encoding="utf-8")
        with pytest.raises(ParseError, match="Empty FASTA header"):
            list(parse_fasta(path))

    def test_raises_for_empty_sequence(self, tmp_path: Path) -> None:
        path = tmp_path / "bad.fasta"
        path.write_text(">seq1\n>seq2\nGGCC\n", encoding="utf-8")
        with pytest.raises(ParseError, match="Empty sequence"):
            list(parse_fasta(path))

    def test_supports_gzip_input(self, gzip_fasta_file: Path) -> None:
        assert len(list(parse_fasta(gzip_fasta_file))) == 2


class TestParseFastq:
    def test_parses_three_records(self, fastq_file: Path) -> None:
        assert len(list(parse_fastq(fastq_file))) == 3

    def test_identifier_is_first_token(self, fastq_file: Path) -> None:
        assert list(parse_fastq(fastq_file))[0].identifier == "read1"

    def test_description_is_rest_of_header(self, fastq_file: Path) -> None:
        assert list(parse_fastq(fastq_file))[0].description == "alpha"

    def test_sequence_is_uppercased(self, tmp_path: Path) -> None:
        path = tmp_path / "lower.fastq"
        path.write_text("@read1\natgc\n+\nIIII\n", encoding="utf-8")
        assert list(parse_fastq(path))[0].sequence == "ATGC"

    def test_max_records_caps_output(self, fastq_file: Path) -> None:
        assert len(list(parse_fastq(fastq_file, max_records=2))) == 2

    def test_blank_lines_between_records_are_skipped(self, tmp_path: Path) -> None:
        path = tmp_path / "blanks.fastq"
        path.write_text("@read1\nATGC\n+\nIIII\n\n@read2\nGGCC\n+\nIIII\n", encoding="utf-8")
        assert len(list(parse_fastq(path))) == 2

    def test_raises_for_missing_file(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            list(parse_fastq(tmp_path / "missing.fastq"))

    def test_raises_for_empty_file(self, empty_file: Path) -> None:
        with pytest.raises(ParseError, match="empty"):
            list(parse_fastq(empty_file))

    def test_raises_for_missing_at_header(self, tmp_path: Path) -> None:
        path = tmp_path / "bad.fastq"
        path.write_text("read1\nATGC\n+\nIIII\n", encoding="utf-8")
        with pytest.raises(ParseError, match="must start with '@'"):
            list(parse_fastq(path))

    def test_raises_for_empty_header(self, tmp_path: Path) -> None:
        path = tmp_path / "bad.fastq"
        path.write_text("@\nATGC\n+\nIIII\n", encoding="utf-8")
        with pytest.raises(ParseError, match="header is empty"):
            list(parse_fastq(path))

    def test_raises_for_missing_separator(self, tmp_path: Path) -> None:
        path = tmp_path / "bad.fastq"
        path.write_text("@read1\nATGC\n-\nIIII\n", encoding="utf-8")
        with pytest.raises(ParseError, match="separator"):
            list(parse_fastq(path))

    def test_raises_for_quality_length_mismatch(self, malformed_fastq_file: Path) -> None:
        with pytest.raises(ParseError, match="quality length"):
            list(parse_fastq(malformed_fastq_file))

    def test_raises_for_incomplete_record(self, incomplete_fastq_file: Path) -> None:
        with pytest.raises(ParseError, match="incomplete"):
            list(parse_fastq(incomplete_fastq_file))

    def test_supports_gzip_input(self, gzip_fastq_file: Path) -> None:
        assert len(list(parse_fastq(gzip_fastq_file))) == 2


class TestParseFile:
    def test_dispatches_fasta(self, fasta_file: Path) -> None:
        assert all(isinstance(record, FastaRecord) for record in parse_file(fasta_file, "fasta"))

    def test_dispatches_fastq(self, fastq_file: Path) -> None:
        assert all(isinstance(record, FastqRecord) for record in parse_file(fastq_file, "fastq"))

    def test_rejects_unknown_file_type(self, fasta_file: Path) -> None:
        with pytest.raises(ValueError, match="Unsupported"):
            list(parse_file(fasta_file, "bam"))
