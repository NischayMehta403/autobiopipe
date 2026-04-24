"""CLI integration tests for the Typer entrypoint."""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from autobiopipe._version import __version__
from autobiopipe.cli import app


runner = CliRunner()


def test_version_command_prints_version() -> None:
    result = runner.invoke(app, ["version"])

    assert result.exit_code == 0
    assert f"AutoBioPipe {__version__}" in result.stdout
    assert "Python " in result.stdout


def test_detect_command_prints_detection_summary(fastq_file: Path) -> None:
    result = runner.invoke(app, ["detect", str(fastq_file)])

    assert result.exit_code == 0
    assert fastq_file.name in result.stdout
    assert "type=fastq" in result.stdout
    assert "sequencing=single-end" in result.stdout


def test_detect_command_returns_exit_code_1_for_missing_file(tmp_path: Path) -> None:
    missing_file = tmp_path / "missing.fastq"

    result = runner.invoke(app, ["detect", str(missing_file)])

    assert result.exit_code == 1
    assert "Detection failed:" in result.stdout


def test_run_command_generates_reports_and_fail_exit_code(fastq_file: Path, tmp_path: Path) -> None:
    output_dir = tmp_path / "cli-run"

    result = runner.invoke(app, ["run", str(fastq_file), "--output-dir", str(output_dir)])

    assert result.exit_code == 2
    assert (output_dir / f"{fastq_file.name}_report.json").exists()
    assert (output_dir / f"{fastq_file.name}_report.csv").exists()
    assert (output_dir / f"{fastq_file.name}_report.pdf").exists()


def test_report_command_renders_saved_json_report(fastq_file: Path, tmp_path: Path) -> None:
    output_dir = tmp_path / "cli-report"
    run_result = runner.invoke(app, ["run", str(fastq_file), "--output-dir", str(output_dir)])

    assert run_result.exit_code == 2

    report_path = output_dir / f"{fastq_file.name}_report.json"
    result = runner.invoke(app, ["report", str(report_path)])

    assert result.exit_code == 0
    assert "Saved AutoBioPipe Report" in result.stdout
    assert "Findings" in result.stdout


def test_report_command_returns_exit_code_1_for_missing_json(tmp_path: Path) -> None:
    missing_report = tmp_path / "missing_report.json"

    result = runner.invoke(app, ["report", str(missing_report)])

    assert result.exit_code == 1
    assert "Report rendering failed:" in result.stdout
