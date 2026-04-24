"""Integration tests for pipeline execution and report generation."""

from __future__ import annotations

import json
from pathlib import Path

from autobiopipe.config import AutoBioPipeConfig, QCConfig
from autobiopipe.pipeline import run_pipeline


def test_pipeline_writes_all_reports(fastq_file: Path, tmp_path: Path) -> None:
    result = run_pipeline(fastq_file, output_dir=tmp_path)
    assert result.output_json.exists()
    assert result.output_csv.exists()
    assert result.output_pdf.exists()
    assert result.visualizations.gc_distribution.exists()
    assert result.visualizations.length_distribution.exists()
    assert result.visualizations.quality_profile.exists()


def test_pipeline_result_contains_stage_timings(fastq_file: Path, tmp_path: Path) -> None:
    result = run_pipeline(fastq_file, output_dir=tmp_path)
    assert {
        "detect",
        "parse",
        "qc",
        "biology",
        "visualization",
        "decision",
        "report",
    } <= set(result.stage_timings)


def test_pipeline_json_report_contains_decision(fastq_file: Path, tmp_path: Path) -> None:
    result = run_pipeline(fastq_file, output_dir=tmp_path)
    payload = json.loads(result.output_json.read_text(encoding="utf-8"))
    assert payload["decision"]["overall_status"] in {"PASS", "WARNING", "FAIL"}
    assert payload["qc"]["biology"]["sequence_scale"]
    assert payload["visualizations"]["gc_distribution"].endswith("gc_distribution.png")


def test_pipeline_fail_status_is_returned(tmp_path: Path) -> None:
    path = tmp_path / "low.fastq"
    path.write_text("@r1\nATGC\n+\n!!!!\n", encoding="utf-8")
    config = AutoBioPipeConfig(qc=QCConfig(min_avg_quality=20.0, min_read_length=1))
    result = run_pipeline(path, output_dir=tmp_path / "out", config=config)
    assert result.decision.overall_status == "FAIL"
