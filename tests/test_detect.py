"""Tests for autobiopipe.detect"""

from pathlib import Path
import pytest

from autobiopipe.detect import DetectionResult, detect_file


class TestDetectFile:
    def test_detects_fasta(self, fasta_file: Path) -> None:
        result = detect_file(fasta_file)
        assert result.file_type == "fasta"
        assert result.confidence == "high"

    def test_detects_fastq(self, fastq_file: Path) -> None:
        result = detect_file(fastq_file)
        assert result.file_type == "fastq"
        assert result.confidence == "high"

    def test_not_found_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            detect_file(tmp_path / "missing.fasta")

    def test_returns_detection_result(self, fasta_file: Path) -> None:
        assert isinstance(detect_file(fasta_file), DetectionResult)

    def test_not_compressed(self, fasta_file: Path) -> None:
        assert detect_file(fasta_file).is_compressed is False

    def test_single_end_default(self, fasta_file: Path) -> None:
        assert detect_file(fasta_file).sequencing_type == "single-end"

    def test_paired_end_detected(self, tmp_path: Path) -> None:
        r1 = tmp_path / "sample_R1.fastq"
        r2 = tmp_path / "sample_R2.fastq"
        content = "@r1\nATGC\n+\nIIII\n"
        r1.write_text(content)
        r2.write_text(content)
        result = detect_file(r1)
        assert result.sequencing_type == "paired-end"
        assert result.paired_mate.name == "sample_R2.fastq"

    def test_unknown_adds_warning(self, tmp_path: Path) -> None:
        p = tmp_path / "weird.xyz"
        p.write_text("not a biological file\n")
        result = detect_file(p)
        assert result.file_type == "unknown"
        assert len(result.warnings) > 0

    def test_summary_contains_format(self, fasta_file: Path) -> None:
        summary = detect_file(fasta_file).summary()
        assert "FASTA" in summary
        assert "single-end" in summary
