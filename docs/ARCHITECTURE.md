# Architecture

## High-Level Flow

AutoBioPipe runs seven stages:

1. `detect`: infer file type, compression, and paired-end hints
2. `parse`: stream FASTA or FASTQ records into typed dataclasses
3. `qc`: compute metrics for GC, lengths, qualities, and N fraction
4. `biology`: infer sequence scale, GC context, coding potential, duplication, and GC skew
5. `visualization`: render scientific plots to PNG files
6. `decision`: apply registry-backed QC rules to the metrics
7. `report`: render terminal, JSON, CSV, and PDF outputs

## Module Relationships

- `config.py` owns threshold and pipeline settings
- `detect.py` provides `DetectionResult`, used by the pipeline and terminal reports
- `parser.py` provides `FastaRecord`, `FastqRecord`, `ParseError`, and streaming parsers
- `qc.py` transforms records into `QCResult`
- `biology.py` enriches `QCResult` with biological interpretation
- `visualization.py` emits plot artifacts used by the reports
- `decision.py` transforms `QCResult` plus config into `DecisionReport`
- `report.py` renders the detection, QC, biology, visualization, and decision outputs
- `pipeline.py` coordinates stage ordering, logging, timings, and report persistence

## Design Choices

- Parsers are generators for memory efficiency.
- QC and report stages use dataclasses for explicit structured outputs.
- Decision rules are registered through a decorator so new rules can be added without editing the engine loop.
- Biological heuristics are intentionally lightweight and interpretable rather than model-driven.
- Visualization is generated with plain matplotlib using a non-interactive backend for CLI and CI safety.
- Config is TOML-based to stay readable for researchers and easy to version-control.
