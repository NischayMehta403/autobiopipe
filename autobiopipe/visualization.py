"""Visualization utilities for AutoBioPipe."""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
import tempfile

from autobiopipe.config import VisualizationConfig
from autobiopipe.parser import AnyRecord, FastqRecord
from autobiopipe.qc import QCResult


@dataclass
class VisualizationArtifacts:
    """Paths to generated visualization files."""

    gc_distribution: Path
    length_distribution: Path
    quality_profile: Path


def _load_pyplot() -> object:
    """Import matplotlib.pyplot with a non-interactive backend."""
    os.environ.setdefault(
        "MPLCONFIGDIR",
        str(Path(tempfile.gettempdir()) / "autobiopipe-mpl-cache"),
    )
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _apply_style(plt: object, config: VisualizationConfig) -> None:
    """Apply a simple scientific plotting style."""
    plt.style.use("default")
    plt.rcParams.update(
        {
            "figure.dpi": config.dpi,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def _save_gc_distribution(
    records: list[AnyRecord],
    output_path: Path,
    config: VisualizationConfig,
) -> Path:
    """Save the GC-distribution histogram."""
    plt = _load_pyplot()
    _apply_style(plt, config)
    values = [record.gc_content() for record in records]
    figure, axis = plt.subplots(figsize=(6.5, 4.0))
    axis.hist(values, bins=min(config.histogram_bins, max(5, len(values))), color="#2A6F97")
    axis.set_title("GC Distribution")
    axis.set_xlabel("GC content (%)")
    axis.set_ylabel("Sequence count")
    figure.tight_layout()
    figure.savefig(output_path)
    plt.close(figure)
    return output_path


def _save_length_distribution(
    records: list[AnyRecord],
    output_path: Path,
    config: VisualizationConfig,
) -> Path:
    """Save the sequence-length histogram."""
    plt = _load_pyplot()
    _apply_style(plt, config)
    values = [len(record) for record in records]
    figure, axis = plt.subplots(figsize=(6.5, 4.0))
    axis.hist(values, bins=min(config.histogram_bins, max(5, len(values))), color="#468FAF")
    axis.set_title("Length Distribution")
    axis.set_xlabel("Length (bp)")
    axis.set_ylabel("Sequence count")
    figure.tight_layout()
    figure.savefig(output_path)
    plt.close(figure)
    return output_path


def _save_quality_profile(
    records: list[AnyRecord],
    qc_result: QCResult,
    output_path: Path,
    config: VisualizationConfig,
) -> Path:
    """Save the quality profile plot or a FASTA placeholder."""
    plt = _load_pyplot()
    _apply_style(plt, config)
    figure, axis = plt.subplots(figsize=(6.5, 4.0))
    if qc_result.quality is not None:
        positions = list(range(1, len(qc_result.quality.per_position_means) + 1))
        axis.plot(positions, qc_result.quality.per_position_means, color="#01497C", linewidth=2.0)
        axis.set_title("Per-base Quality Profile")
        axis.set_xlabel("Base position")
        axis.set_ylabel("Mean Phred score")
    else:
        axis.text(
            0.5,
            0.5,
            "Quality profile not available for FASTA input",
            ha="center",
            va="center",
            fontsize=11,
            transform=axis.transAxes,
        )
        axis.set_title("Per-base Quality Profile")
        axis.set_xticks([])
        axis.set_yticks([])
    figure.tight_layout()
    figure.savefig(output_path)
    plt.close(figure)
    return output_path


def generate_visualizations(
    records: list[AnyRecord],
    qc_result: QCResult,
    output_dir: Path,
    config: VisualizationConfig,
) -> VisualizationArtifacts:
    """Generate all required AutoBioPipe figures.

    Args:
        records: Parsed FASTA or FASTQ records.
        qc_result: Computed QC metrics.
        output_dir: Directory to write plots into.
        config: Visualization configuration.

    Returns:
        VisualizationArtifacts containing all output paths.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    gc_path = _save_gc_distribution(records, output_dir / "gc_distribution.png", config)
    length_path = _save_length_distribution(
        records,
        output_dir / "length_distribution.png",
        config,
    )
    quality_path = _save_quality_profile(
        records,
        qc_result,
        output_dir / "quality_profile.png",
        config,
    )
    return VisualizationArtifacts(
        gc_distribution=gc_path,
        length_distribution=length_path,
        quality_profile=quality_path,
    )
