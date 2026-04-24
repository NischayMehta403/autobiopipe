"""Tests for input-file detection."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from autobiopipe.detect import DetectionResult, detect_file


class TestDetectFile:
    def test_detects_fasta_by_content(self, fasta_file: Path) -> None:
        result = detect_file(fasta_file)
        assert result.file_type == "fasta"
        assert result.confidence == "high"

    def test_detects_fastq_by_content(self, fastq_file: Path) -> None:
        result = detect_file(fastq_file)
        assert result.file_type == "fastq"
        assert result.confidence == "high"

    def test_detects_gzip_fastq(self, gzip_fastq_file: Path) -> None:
        result = detect_file(gzip_fastq_file)
        assert result.file_type == "fastq"
        assert result.is_compressed is True

    def test_falls_back_to_extension_for_unknown_content(self, tmp_path: Path) -> None:
        path = tmp_path / "reads.fastq"
        path.write_text("not_a_real_record\n", encoding="utf-8")
        result = detect_file(path)
        assert result.file_type == "fastq"
        assert result.confidence == "medium"

    def test_warns_when_extension_disagrees_with_content(self, tmp_path: Path) -> None:
        path = tmp_path / "reads.fasta"
        path.write_text("@read1\nATGC\n+\nIIII\n", encoding="utf-8")
        result = detect_file(path)
        assert result.file_type == "fastq"
        assert result.warnings

    def test_detects_single_end_by_default(self, fasta_file: Path) -> None:
        assert detect_file(fasta_file).sequencing_type == "single-end"

    def test_detects_paired_end_when_mate_exists(self, paired_fastq_files: tuple[Path, Path]) -> None:
        r1, r2 = paired_fastq_files
        result = detect_file(r1)
        assert result.sequencing_type == "paired-end"
        assert result.paired_mate == r2.resolve()

    def test_warns_when_paired_mate_is_missing(self, tmp_path: Path) -> None:
        path = tmp_path / "sample_R1.fastq"
        path.write_text("@read1\nATGC\n+\nIIII\n", encoding="utf-8")
        result = detect_file(path)
        assert result.sequencing_type == "paired-end"
        assert result.paired_mate is None
        assert result.warnings

    def test_returns_detection_result_dataclass(self, fasta_file: Path) -> None:
        assert isinstance(detect_file(fasta_file), DetectionResult)

    def test_reports_unknown_for_empty_file(self, empty_file: Path) -> None:
        result = detect_file(empty_file)
        assert result.file_type == "fasta"

    def test_does_not_treat_accession_version_as_paired_end(self, tmp_path: Path) -> None:
        path = tmp_path / "genbank_J01673.1.fasta"
        path.write_text(">J01673.1 example\nATGC\n", encoding="utf-8")
        result = detect_file(path)
        assert result.file_type == "fasta"
        assert result.sequencing_type == "single-end"
        assert result.paired_mate is None
        assert not result.warnings

    def test_unknown_content_and_unknown_extension(self, tmp_path: Path) -> None:
        path = tmp_path / "reads.bin"
        path.write_text("???\n", encoding="utf-8")
        result = detect_file(path)
        assert result.file_type == "unknown"
        assert result.warnings

    def test_summary_mentions_confidence(self, fastq_file: Path) -> None:
        assert "confidence" in detect_file(fastq_file).summary()
