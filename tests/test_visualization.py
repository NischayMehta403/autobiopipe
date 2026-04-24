"""Tests for visualization output generation."""

from __future__ import annotations

from pathlib import Path

from autobiopipe.config import VisualizationConfig
from autobiopipe.parser import parse_fastq
from autobiopipe.qc import run_qc
from autobiopipe.visualization import generate_visualizations


def test_generate_visualizations_creates_png_files(
    fastq_file: Path,
    tmp_path: Path,
) -> None:
    records = list(parse_fastq(fastq_file))
    qc_result = run_qc(records, fastq_file.name, "fastq", quality_threshold=20.0)
    artifacts = generate_visualizations(
        records,
        qc_result,
        tmp_path,
        VisualizationConfig(enabled=True, dpi=100, histogram_bins=10),
    )
    assert artifacts.gc_distribution.exists()
    assert artifacts.length_distribution.exists()
    assert artifacts.quality_profile.exists()
