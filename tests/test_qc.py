"""Tests for autobiopipe.qc"""

from pathlib import Path
from typing import List
import pytest

from autobiopipe.parser import FastaRecord, FastqRecord
from autobiopipe.qc import (
    compute_gc_stats,
    compute_length_stats,
    compute_quality_stats,
    compute_n_fraction,
    run_qc,
)


def _fasta(seqs: List[str]) -> List[FastaRecord]:
    return [FastaRecord(f"s{i}", f"s{i}", s) for i, s in enumerate(seqs)]


def _fastq(seqs: List[str], quals: List[str]) -> List[FastqRecord]:
    return [FastqRecord(f"r{i}", f"r{i}", s, q)
            for i, (s, q) in enumerate(zip(seqs, quals))]


# ── GC stats ──────────────────────────────────────────────────


class TestGCStats:
    def test_pure_gc(self) -> None:
        stats = compute_gc_stats(_fasta(["GCGC", "GCGC"]))
        assert stats.mean == pytest.approx(100.0)

    def test_balanced(self) -> None:
        stats = compute_gc_stats(_fasta(["ATGC"]))
        assert stats.mean == pytest.approx(50.0)

    def test_empty(self) -> None:
        stats = compute_gc_stats([])
        assert stats.mean == 0.0

    def test_is_normal_true(self) -> None:
        stats = compute_gc_stats(_fasta(["ATGCATGC"]))
        assert stats.is_normal(35.0, 65.0) is True

    def test_is_normal_false(self) -> None:
        stats = compute_gc_stats(_fasta(["AAAAAAAAAAAAA"]))
        assert stats.is_normal(35.0, 65.0) is False

    def test_histogram_sums_to_record_count(self) -> None:
        records = _fasta(["ATGC", "GCGC", "ATAT"])
        stats = compute_gc_stats(records)
        assert sum(stats.histogram.values()) == 3


# ── Length stats ──────────────────────────────────────────────


class TestLengthStats:
    def test_uniform_lengths(self) -> None:
        stats = compute_length_stats(_fasta(["ATGC", "GCTA"]))
        assert stats.mean == pytest.approx(4.0)
        assert stats.coefficient_of_variation == pytest.approx(0.0)

    def test_total_bases(self) -> None:
        stats = compute_length_stats(_fasta(["ATGC", "GCGC"]))
        assert stats.total_bases == 8

    def test_min_max(self) -> None:
        stats = compute_length_stats(_fasta(["AT", "ATGCATGC"]))
        assert stats.min_value == 2
        assert stats.max_value == 8

    def test_empty(self) -> None:
        stats = compute_length_stats([])
        assert stats.total_bases == 0


# ── Quality stats ─────────────────────────────────────────────


class TestQualityStats:
    def test_high_quality(self) -> None:
        records = _fastq(["ATGC", "GCTA"], ["IIII", "IIII"])
        stats = compute_quality_stats(records, quality_threshold=20.0)
        assert stats.mean_per_read == pytest.approx(40.0)
        assert stats.reads_below_threshold == 0

    def test_low_quality(self) -> None:
        records = _fastq(["ATGC", "GCTA"], ["!!!!", "!!!!"])
        stats = compute_quality_stats(records, quality_threshold=20.0)
        assert stats.reads_below_threshold == 2
        assert stats.reads_below_threshold_pct == pytest.approx(100.0)

    def test_per_position_length(self) -> None:
        records = _fastq(["ATGC"], ["IIII"])
        stats = compute_quality_stats(records)
        assert len(stats.per_position_mean) == 4

    def test_empty(self) -> None:
        stats = compute_quality_stats([])
        assert stats.mean_per_read == 0.0


# ── N fraction ────────────────────────────────────────────────


class TestNFraction:
    def test_no_ns(self) -> None:
        assert compute_n_fraction(_fasta(["ATGC"])) == pytest.approx(0.0)

    def test_all_ns(self) -> None:
        assert compute_n_fraction(_fasta(["NNNN"])) == pytest.approx(1.0)

    def test_half_ns(self) -> None:
        assert compute_n_fraction(_fasta(["ATNN"])) == pytest.approx(0.5)

    def test_empty(self) -> None:
        assert compute_n_fraction([]) == pytest.approx(0.0)


# ── run_qc integration ────────────────────────────────────────


class TestRunQC:
    def test_fasta_has_no_quality(self, fasta_file: Path) -> None:
        from autobiopipe.parser import parse_fasta
        records = list(parse_fasta(fasta_file))
        result = run_qc(records, "sample.fasta", "fasta")
        assert result.quality is None

    def test_fastq_has_quality(self, fastq_file: Path) -> None:
        from autobiopipe.parser import parse_fastq
        records = list(parse_fastq(fastq_file))
        result = run_qc(records, "sample.fastq", "fastq")
        assert result.quality is not None

    def test_high_n_adds_warning(self, tmp_path: Path) -> None:
        from autobiopipe.parser import parse_fasta
        p = tmp_path / "n.fasta"
        p.write_text(">s1\n" + "N" * 50 + "A" * 5 + "\n")
        records = list(parse_fasta(p))
        result = run_qc(records, "n.fasta", "fasta")
        assert any("N-base" in w for w in result.warnings)
