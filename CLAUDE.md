# CLAUDE.md

## Purpose

This repository implements AutoBioPipe, a modular CLI for FASTA/FASTQ quality control. The primary user flow is:

1. Detect input characteristics in `autobiopipe.detect`
2. Parse records in `autobiopipe.parser`
3. Compute metrics in `autobiopipe.qc`
4. Infer biological context in `autobiopipe.biology`
5. Generate figures in `autobiopipe.visualization`
6. Apply rules in `autobiopipe.decision`
7. Emit reports in `autobiopipe.report`

## Architecture

- `autobiopipe.config`: Typed config dataclasses plus TOML loading and validation
- `autobiopipe.detect`: Content sniffing, extension fallback, gzip handling, paired-end heuristics
- `autobiopipe.parser`: Streaming FASTA/FASTQ parsing with `ParseError`
- `autobiopipe.qc`: Aggregated metrics dataclasses and QC orchestration
- `autobiopipe.biology`: GC interpretation, sequence-scale inference, ORF heuristics, coding-potential scoring, complexity metrics
- `autobiopipe.visualization`: Matplotlib-based plot generation for GC, length, and quality summaries
- `autobiopipe.decision`: Registry-based rules `QC001` to `QC015` with weighted scoring
- `autobiopipe.pipeline`: End-to-end orchestration with logging and timings
- `autobiopipe.report`: Terminal, JSON, CSV, and PDF outputs
- `autobiopipe.ai_explain`: Placeholder for optional Gemini integration
- `autobiopipe.ml`: Reserved for future ML-driven classifiers

## Extension Patterns

- Add a new decision rule by writing a function in `autobiopipe.decision` with signature `(QCResult, AutoBioPipeConfig) -> list[Finding]` and decorating it with `@register_rule("QCXYZ")`.
- Extend config by adding fields to the config dataclasses and merging them in `load_config`.
- Keep parsers generator-based. Avoid reading large files into memory unless a downstream stage explicitly needs materialization.

## Gotchas

- `report.py` should stay JSON-safe. Convert dataclasses, enums, and paths before serialization.
- `ParseError` should include 1-based line numbers whenever possible.
- FASTQ quality scores are ASCII-33 encoded.
- Paired-end detection is heuristic and filename-based. Do not treat it as authoritative metadata.
- `visualization.py` should keep matplotlib imports lazy and use a writable `MPLCONFIGDIR` so CLI runs work in constrained environments.
- `CODEX.md` must be updated whenever behavior, code, tests, docs, or verification changes.

## Expected Quality Bar

- Google-style docstrings on public functions
- Type hints on all function signatures
- Config-driven thresholds in the decision engine
- Tests for both success and failure paths
- Avoid hidden side effects and keep CLI behavior deterministic
