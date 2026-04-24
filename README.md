# AutoBioPipe

[![Tests](https://github.com/NischayMehta403/autobiopipe/actions/workflows/test.yml/badge.svg)](https://github.com/NischayMehta403/autobiopipe/actions/workflows/test.yml)
[![Python](https://img.shields.io/badge/python-3.9%2B-1f6feb.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-0f766e.svg)](LICENSE)

AutoBioPipe is a Python CLI for sequence-file quality control that goes beyond basic FASTA/FASTQ validation.

It detects file structure, computes core QC metrics, adds lightweight biological interpretation, renders plots, and produces terminal, JSON, CSV, and PDF reports in one run.

## Why AutoBioPipe

Most quick QC scripts stop at counts and averages. AutoBioPipe is designed to give you a fuller first-pass read on a dataset:

- format detection for FASTA, FASTQ, and `.gz` inputs
- streaming parsers with explicit error handling
- GC, length, quality, and N-content QC metrics
- lightweight biological heuristics for coding potential, GC context, duplication, low complexity, homopolymers, and GC skew
- generated plots for GC distribution, read lengths, and quality profiles
- a configurable decision engine with rule-based findings and severity scoring
- machine-readable and human-readable reports from the same run

## Pipeline At A Glance

AutoBioPipe runs the following stages:

`detect -> parse -> qc -> biology -> visualization -> decision -> report`

That means a single command can tell you:

- what kind of file you gave it
- whether the records parse cleanly
- whether basic QC metrics look suspicious
- whether the sequences show biologically unusual patterns
- what concrete findings were triggered

## Installation

```bash
git clone https://github.com/NischayMehta403/autobiopipe.git
cd autobiopipe
pip install -e .
```

For development:

```bash
pip install -e .[dev]
```

## Quick Start

Run the full pipeline on the bundled sample FASTQ:

```bash
autobio run examples/sample.fastq --output-dir results
```

Inspect detection only:

```bash
autobio detect examples/sample.fasta
```

Render a saved JSON report back to the terminal:

```bash
autobio report results/sample.fastq_report.json
```

Show the installed version:

```bash
autobio version
```

## What You Get

A typical run writes:

- `*_report.json`
- `*_report.csv`
- `*_report.pdf`
- `gc_distribution.png`
- `length_distribution.png`
- `quality_profile.png`

The JSON report is the canonical structured output. CSV is useful for spreadsheets and downstream tabular review. PDF is intended for quick sharing and archival summaries.

## CLI

`autobio run <input_file> [--output-dir DIR] [--config FILE] [--verbose]`

- runs the complete pipeline
- returns exit code `0` for `PASS` or `WARNING`
- returns exit code `1` for runtime or input errors
- returns exit code `2` when the pipeline completes but the overall decision is `FAIL`

`autobio detect <input_file>`

- prints detected file type, compression status, confidence, and inferred sequencing layout

`autobio report <json_file>`

- re-renders a previously generated JSON report in the terminal

`autobio version`

- prints the AutoBioPipe version and Python version

## Configuration

AutoBioPipe uses TOML configuration files.

Examples:

- [examples/config.toml](examples/config.toml)
- [examples/config_strict.toml](examples/config_strict.toml)

Main config sections:

- `[qc]` for core QC thresholds
- `[biology]` for biological heuristics and advanced rule thresholds
- `[visualization]` for plot settings
- `[pipeline]` for runtime and output behavior

Run with a custom config:

```bash
autobio run examples/sample.fastq --config examples/config.toml
```

## Project Structure

- [autobiopipe/](autobiopipe) contains the CLI and pipeline modules
- [tests/](tests) contains parser, QC, decision, pipeline, visualization, biology, and CLI coverage
- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) explains the stage design
- [docs/CONTRIBUTING.md](docs/CONTRIBUTING.md) documents contribution guidance

## Verification

The test suite is exercised with `pytest`, and the repository includes GitHub Actions CI in [.github/workflows/test.yml](.github/workflows/test.yml).

## Notes

- Empty files, malformed FASTQ records, multiline FASTA sequences, gzip compression, and paired-end filename hints are handled explicitly.
- The bundled example files are intentionally small so the project stays lightweight for demos and tests.
- `examples/sample.fastq` is intentionally harsh under the advanced rule set and may return `FAIL`; treat it as a regression fixture, not as a biologically representative "good" sample.
