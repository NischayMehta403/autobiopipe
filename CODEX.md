# CODEX.md

## Purpose

This file is a session handoff for Codex so a new chat can resume work on AutoBioPipe without rebuilding context from scratch.

## Session Preference

User instruction to preserve across future work:

- every time code, tests, docs, config, behavior, or verification state changes, update `CODEX.md` in the same pass

## Repository Location

- Root repo: `/home/Fred/autobiopipe`
- Main package: `/home/Fred/autobiopipe/autobiopipe`

## Current Project State

AutoBioPipe has been scaffolded and implemented as a production-ready Python CLI project for FASTA/FASTQ QC.

Implemented files and areas:

- Packaging:
  - `pyproject.toml`
  - `requirements.txt`
  - `.gitignore`
  - `LICENSE`
- Docs:
  - `README.md`
  - `CLAUDE.md`
  - `docs/CONTRIBUTING.md`
  - `docs/ARCHITECTURE.md`
- Examples:
  - `examples/sample.fasta`
  - `examples/sample.fastq`
  - `examples/config.toml`
  - `examples/config_strict.toml`
- CI:
  - `.github/workflows/test.yml`
- Package:
  - `autobiopipe/__init__.py`
  - `autobiopipe/_version.py`
  - `autobiopipe/cli.py`
  - `autobiopipe/biology.py`
  - `autobiopipe/config.py`
  - `autobiopipe/detect.py`
  - `autobiopipe/parser.py`
  - `autobiopipe/qc.py`
  - `autobiopipe/decision.py`
  - `autobiopipe/pipeline.py`
  - `autobiopipe/report.py`
  - `autobiopipe/visualization.py`
  - `autobiopipe/ai_explain.py`
  - `autobiopipe/ml/__init__.py`
- Tests:
  - `tests/test_biology.py`
  - `tests/conftest.py`
  - `tests/test_parser.py`
  - `tests/test_qc.py`
  - `tests/test_detect.py`
  - `tests/test_decision.py`
  - `tests/test_pipeline.py`
  - `tests/test_visualization.py`

## Implemented Behavior

### CLI

Commands implemented in `autobiopipe/cli.py`:

- `autobio run`
- `autobio detect`
- `autobio report`
- `autobio version`

Exit behavior:

- `0`: success / pass / warning
- `1`: runtime or input error
- `2`: pipeline completed but overall status is `FAIL`

### Config

Config loader implemented in `autobiopipe/config.py`.

Supported config sections:

- `[qc]`
- `[biology]`
- `[visualization]`
- `[pipeline]`

Key QC thresholds:

- `min_avg_quality`
- `gc_content_min`
- `gc_content_max`
- `min_read_length`
- `length_cv_threshold`
- `n_fraction_warning_threshold`
- `n_fraction_critical_threshold`
- `small_dataset_threshold`
- `enable_small_dataset_warning`

Pipeline settings:

- `verbose`
- `output_dir`
- `max_records`

Biology settings include thresholds for:

- sequence scale inference
- coding potential
- low complexity
- homopolymers
- duplication
- GC extremes
- GC skew
- short sequence fractions

Visualization settings include:

- `enabled`
- `dpi`
- `histogram_bins`

### Detection

Implemented in `autobiopipe/detect.py`:

- content sniffing using first non-empty character
- extension fallback
- gzip detection
- paired-end filename heuristics
- confidence and warning reporting

### Parser

Implemented in `autobiopipe/parser.py`:

- `FastaRecord`
- `FastqRecord`
- `ParseError`
- `parse_fasta()`
- `parse_fastq()`
- `parse_file()`

Important behavior:

- generator-based parsing
- gzip-transparent input
- FASTQ ASCII-33 quality decoding
- line-numbered parse errors where applicable
- multiline FASTA support

### QC

Implemented in `autobiopipe/qc.py`:

- `compute_gc_stats()`
- `compute_length_stats()`
- `compute_quality_stats()`
- `compute_n_fraction()`
- `run_qc()`

Metrics include warnings for anomalies.

`QCResult` now also carries:

- `sequence_scale`
- `biology`

where `biology` is a structured biological summary produced in the new biology stage.

### Biology

Implemented in `autobiopipe/biology.py`:

- `GCReferenceDatabase`
- `infer_organism_type()`
- `evaluate_gc_context()`
- `classify_sequence_scale()`
- `estimate_coding_potential()`
- `analyze_biology()`

Biological summary currently includes:

- inferred organism/context type
- GC context and explanation
- sequence scale
- coding potential classification and score
- ORF presence fraction
- codon-usage bias estimate
- low-complexity fraction
- homopolymer fraction
- duplication fraction
- GC skew mean absolute value
- short-sequence fraction

### Visualization

Implemented in `autobiopipe/visualization.py`:

- `generate_visualizations()`
- `VisualizationArtifacts`

Generated files:

- `gc_distribution.png`
- `length_distribution.png`
- `quality_profile.png`

### Decision Engine

Implemented in `autobiopipe/decision.py`.

Built-in rules:

- `QC001`: low quality
- `QC002`: GC anomaly
- `QC003`: length consistency / minimum length
- `QC004`: N fraction
- `QC005`: small dataset
- `QC006`: extreme GC content
- `QC007`: low complexity sequence
- `QC008`: homopolymer runs
- `QC009`: abnormal length distribution
- `QC010`: high variance in read length
- `QC011`: coding potential anomaly
- `QC012`: ORF absence warning
- `QC013`: suspiciously short sequences
- `QC014`: sequence duplication detection
- `QC015`: GC skew imbalance

Pattern used:

- registry via `@register_rule("QCXYZ")`
- `DecisionEngine.evaluate(qc_result)`
- weighted severity scoring:
  - `INFO = 1`
  - `WARNING = 2`
  - `CRITICAL = 3`

### Reporting

Implemented in `autobiopipe/report.py`:

- terminal reporting with Rich
- JSON reporting
- CSV reporting
- PDF reporting
- biological interpretation section
- plot references in all report formats

Important note:

- A minimal fallback PDF writer was added so pipeline runs still produce a PDF even if `reportlab` is not installed in the current environment.
- matplotlib is configured lazily with a writable temp cache directory to avoid runtime noise in constrained environments.

## Test Status

Test suite currently passes.

Last verified command:

```bash
pytest -o addopts=''
```

Observed result:

- `106 passed`

Note:

- `pytest` without overriding `addopts` failed locally because `pytest-cov` was not installed in the current system Python environment, although coverage config remains in `pyproject.toml` for proper dev/CI environments.

## CLI Verification

Verified using the local virtualenv:

```bash
.venv/bin/python -m autobiopipe.cli version
.venv/bin/python -m autobiopipe.cli detect examples/sample.fastq
.venv/bin/python -m autobiopipe.cli detect real_data/genbank_J01673.1.fasta
.venv/bin/python -m autobiopipe.cli run examples/sample.fastq --output-dir tmp_portfolio_run2
```

These worked.

Notes from CLI verification:

- the visualization layer generated PNG plots successfully
- matplotlib cache warnings were fixed by forcing a writable temp `MPLCONFIGDIR`
- `examples/sample.fastq` now returns exit code `2` because the advanced rule set intentionally flags the tiny demo dataset as biologically suspicious

System Python did not have all dependencies installed, so `python -m autobiopipe.cli ...` failed outside the virtualenv due to missing `typer`.

## Real Data Run

A real NCBI GenBank FASTA record was downloaded and analyzed.

Source used:

- NCBI E-utilities `efetch`
- accession: `J01673.1`

Downloaded file:

- `real_data/genbank_J01673.1.fasta`

Generated reports:

- `real_data/results/genbank_J01673.1.fasta_report.json`
- `real_data/results/genbank_J01673.1.fasta_report.csv`
- `real_data/results/genbank_J01673.1.fasta_report.pdf`

Run command:

```bash
.venv/bin/python -m autobiopipe.cli run real_data/genbank_J01673.1.fasta --output-dir real_data/results
```

Observed QC summary:

- file type: `fasta`
- total records: `1`
- length mean: `1880`
- GC mean: `48.56%`
- N fraction: `0.00%`
- overall status: `WARNING`

Triggered finding:

- `QC005` because the dataset had only one record

Latest rerun after the paired-end fix:

- `.venv/bin/python -m autobiopipe.cli run real_data/genbank_J01673.1.fasta --output-dir real_data/results`
- regenerated JSON, CSV, and PDF reports
- detection in saved reports now shows:
  - `sequencing_type = single-end`
  - `warnings = []`

Verified report outputs:

- `real_data/results/genbank_J01673.1.fasta_report.json`
- `real_data/results/genbank_J01673.1.fasta_report.csv`
- `real_data/results/genbank_J01673.1.fasta_report.pdf`

## Fixed Since Initial Handoff

### Paired-end detection false positive on accession suffixes

This bug has been fixed.

Previous issue:

- `autobiopipe/detect.py` interpreted accession version suffixes like `.1` / `.2` as read-pair markers.
- `genbank_J01673.1.fasta` was incorrectly labeled `paired-end`.

Fix applied:

- paired-end filename matching was tightened so:
  - `R1` / `R2` tokens still match with `.`, `_`, or `-`
  - plain `1` / `2` tokens only match with `_` or `-`
- accession suffixes like `.1` no longer trigger paired-end detection

Regression coverage added:

- `genbank_J01673.1.fasta` now tests as `single-end`
- existing `sample_R1.fastq` / `sample_R2.fastq` paired-end behavior still passes

Current real-data behavior:

- `.venv/bin/python -m autobiopipe.cli detect real_data/genbank_J01673.1.fasta`
- now reports `sequencing=single-end`

Pipeline verification after the fix:

- `pytest -o addopts='' tests/test_detect.py tests/test_pipeline.py`
- result: `17 passed`

## Latest Extension Pass

New capabilities added in the latest session:

- biological intelligence layer via `autobiopipe/biology.py`
- visualization layer via `autobiopipe/visualization.py`
- `matplotlib` dependency added to packaging
- advanced decision engine expanded from 5 rules to 15 rules
- weighted severity scoring in the decision report
- biological interpretation and plot references integrated into reports
- pipeline stages expanded to:
  - detect
  - parse
  - qc
  - biology
  - visualization
  - decision
  - report

Latest verification commands:

```bash
.venv/bin/python -m pytest -o addopts='' tests/test_biology.py tests/test_visualization.py tests/test_decision.py tests/test_pipeline.py
.venv/bin/python -m pytest -o addopts=''
.venv/bin/python -m autobiopipe.cli run examples/sample.fastq --output-dir tmp_portfolio_run2
```

Observed results:

- targeted extension tests: `23 passed`
- full suite: `106 passed`
- end-to-end CLI run: successful report and plot generation, with expected exit code `2` on the tiny sample FASTQ because the richer rule set marks it as biologically suspect

## Git / Repo Notes

Important cleanup already performed:

- the nested `autobiopipe/` gitlink/submodule-style entry was converted into normal tracked files in the root repository
- tracked `tests/__pycache__` bytecode files were removed from git tracking

Current repository should now behave like a normal single repository.

Repository naming note:

- root project folder name `autobiopipe/` and Python package name `autobiopipe/` are the same
- this is normal for Python projects and is not a GitHub problem
- current layout is acceptable for pushing to GitHub

## Suggested Next Tasks

If resuming in a new session, good next steps are:

1. Optionally run AutoBioPipe on a larger real multi-record FASTA or a real FASTQ dataset.
2. Optionally recalibrate thresholds or example data so the bundled demo FASTQ is less harsh under the advanced rule set.
3. Optionally verify `pip install -e .` in a clean environment.
4. Optionally add CLI tests if desired.

## Publish Status

User instruction for publish scope:

- include the full current worktree when pushing to GitHub, including generated outputs and supporting artifacts currently present in the repository tree

Current publish action:

- preparing to stage, commit, and push the full repository state to `origin/main`

## How To Resume

In a new session, provide this file plus the repo and say something like:

> Resume from `CODEX.md` and continue from the next suggested task. Keep `CODEX.md` updated with every change.

That should be enough to continue without reconstructing the prior conversation.
