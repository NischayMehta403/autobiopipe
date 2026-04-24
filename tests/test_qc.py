"""Tests for QC metric calculations."""

from __future__ import annotations

from pathlib import Path

import pytest

from autobiopipe.parser import FastaRecord, FastqRecord, parse_fasta, parse_fastq
from autobiopipe.qc import (
    compute_gc_stats,
    compute_length_stats,
    compute_n_fraction,
    compute_quality_stats,
    run_qc,
)


def _fasta(*sequences: str) -> list[FastaRecord]:
    return [FastaRecord(f"seq{i}", "", sequence) for i, sequence in enumerate(sequences, start=1)]


def _fastq(sequences: list[str], qualities: list[str]) -> list[FastqRecord]:
    return [
        FastqRecord(f"read{i}", "", sequence, quality)
        for i, (sequence, quality) in enumerate(zip(sequences, qualities), start=1)
    ]


class TestGCStats:
    def test_empty_records_return_zeroed_stats(self) -> None:
        stats = compute_gc_stats([])
        assert stats.mean == pytest.approx(0.0)
        assert stats.histogram == {}

    def test_balanced_gc_mean(self) -> None:
        assert compute_gc_stats(_fasta("ATGC")).mean == pytest.approx(50.0)

    def test_extreme_gc_warning(self) -> None:
        assert compute_gc_stats(_fasta("GGGG", "GGGG")).warnings

    def test_histogram_counts_all_records(self) -> None:
        stats = compute_gc_stats(_fasta("ATGC", "GGGG", "AAAA"))
        assert sum(stats.histogram.values()) == 3

    def test_min_and_max_values(self) -> None:
        stats = compute_gc_stats(_fasta("AAAA", "GGGG"))
        assert stats.min_value == pytest.approx(0.0)
        assert stats.max_value == pytest.approx(100.0)

    def test_std_dev_nonzero_when_values_vary(self) -> None:
        assert compute_gc_stats(_fasta("AAAA", "GGGG")).std_dev > 0


class TestLengthStats:
    def test_empty_records_return_zeroed_stats(self) -> None:
        stats = compute_length_stats([])
        assert stats.total_bases == 0
        assert stats.warnings

    def test_total_bases(self) -> None:
        assert compute_length_stats(_fasta("ATGC", "GG")).total_bases == 6

    def test_mean_length(self) -> None:
        assert compute_length_stats(_fasta("ATGC", "GGCC")).mean == pytest.approx(4.0)

    def test_median_length(self) -> None:
        assert compute_length_stats(_fasta("A", "AT", "ATG")).median == pytest.approx(2.0)

    def test_cv_zero_for_uniform_lengths(self) -> None:
        assert compute_length_stats(_fasta("ATGC", "GGCC")).coefficient_of_variation == pytest.approx(0.0)

    def test_short_record_warning(self) -> None:
        assert compute_length_stats(_fasta("ATGC")).warnings

    def test_variable_length_warning(self) -> None:
        assert compute_length_stats(_fasta("A", "ATGCATGCATGC")).warnings


class TestQualityStats:
    def test_empty_records_return_zeroed_stats(self) -> None:
        stats = compute_quality_stats([])
        assert stats.mean_per_read == pytest.approx(0.0)
        assert stats.warnings

    def test_mean_per_read(self) -> None:
        stats = compute_quality_stats(_fastq(["ATGC"], ["IIII"]))
        assert stats.mean_per_read == pytest.approx(40.0)

    def test_reads_below_threshold(self) -> None:
        stats = compute_quality_stats(_fastq(["ATGC", "ATGC"], ["!!!!", "IIII"]), quality_threshold=20.0)
        assert stats.reads_below_threshold == 1

    def test_reads_below_threshold_pct(self) -> None:
        stats = compute_quality_stats(_fastq(["ATGC", "ATGC"], ["!!!!", "IIII"]), quality_threshold=20.0)
        assert stats.reads_below_threshold_pct == pytest.approx(50.0)

    def test_per_position_means(self) -> None:
        stats = compute_quality_stats(_fastq(["ATGC"], ["IIII"]))
        assert stats.per_position_means == [40.0, 40.0, 40.0, 40.0]

    def test_quality_warning_when_mean_below_threshold(self) -> None:
        assert compute_quality_stats(_fastq(["ATGC"], ["!!!!"]), quality_threshold=20.0).warnings

    def test_quality_warning_when_many_reads_below_threshold(self) -> None:
        stats = compute_quality_stats(_fastq(["ATGC"] * 10, ["!!!!"] * 2 + ["IIII"] * 8), quality_threshold=20.0)
        assert any("10%" in warning for warning in stats.warnings)


class TestNFraction:
    def test_empty_records_have_zero_fraction(self) -> None:
        assert compute_n_fraction([]) == pytest.approx(0.0)

    def test_no_ns(self) -> None:
        assert compute_n_fraction(_fasta("ATGC")) == pytest.approx(0.0)

    def test_all_ns(self) -> None:
        assert compute_n_fraction(_fasta("NNNN")) == pytest.approx(1.0)

    def test_half_ns(self) -> None:
        assert compute_n_fraction(_fasta("ATNN")) == pytest.approx(0.5)


class TestRunQC:
    def test_run_qc_for_fasta_has_no_quality(self, fasta_file: Path) -> None:
        result = run_qc(list(parse_fasta(fasta_file)), fasta_file.name, "fasta")
        assert result.quality is None

    def test_run_qc_for_fastq_has_quality(self, fastq_file: Path) -> None:
        result = run_qc(list(parse_fastq(fastq_file)), fastq_file.name, "fastq")
        assert result.quality is not None

    def test_total_records_matches_input(self, fasta_file: Path) -> None:
        result = run_qc(list(parse_fasta(fasta_file)), fasta_file.name, "fasta")
        assert result.total_records == 3

    def test_high_n_fraction_adds_warning(self) -> None:
        result = run_qc(_fasta("NNNN", "AANN"), "n.fasta", "fasta")
        assert any("Ambiguous base fraction" in warning for warning in result.warnings)

    def test_empty_records_add_warning(self) -> None:
        result = run_qc([], "empty.fastq", "fastq")
        assert any("no parseable records" in warning.lower() for warning in result.warnings)

    def test_quality_threshold_is_propagated(self) -> None:
        result = run_qc(_fastq(["ATGC"], ["5555"]), "reads.fastq", "fastq", quality_threshold=25.0)
        assert result.quality is not None
        assert result.quality.threshold_used == pytest.approx(25.0)
