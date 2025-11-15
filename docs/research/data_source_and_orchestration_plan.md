# ABCA4 Campaign – Data & Package Verification Plan

Purpose: hand this to another agent so they can (1) confirm where every required dataset lives, (2) capture the latest safe package versions we should target, and (3) outline how all moving pieces get orchestrated with those packages.

---

## 1. Objectives & Deliverables

1. **Data audit**: Verified locations, access paths, and download recipes for every dataset we need (ClinVar VCF, gnomAD subset, DECIPHER optional, SpliceAI/AlphaMissense tables, domain annotations, DNA FM weights).
2. **Package matrix**: Table of up-to-date versions + install commands for genomics IO, model inference, Strand SDK deps, and orchestration tooling, with verification dates.
3. **Orchestration brief**: Recommendation on how to stitch data prep → feature computation → Strand optimization using the chosen packages (CLI scripts, notebooks, pipelines, tracking).

Output format: a Markdown report dropped back into this folder (`abca4_campaign/`) plus any helper scripts/notebooks under `scripts/`.

---

## 2. Data Verification Tasks

1. **ClinVar VUS backlog**
   - Confirm latest GRCh38 VCF endpoint (FTP + release date).
   - Document `wget`/`curl` command, expected checksum, and gzip size.
   - Define filtering schema (GENEINFO contains ABCA4 + CLNSIG ∈ {Uncertain, Conflicting}).

2. **gnomAD allele frequencies**
   - Locate per-variant TSV or VCF slices for ABCA4.
   - Capture whether to use gnomAD v4.1.0 (preferred) and note liftover needs.
   - Specify method (bcftools view vs gnomAD API) and columns to retain.

3. **DECIPHER optional feed**
   - Check data use terms + API endpoints for ABCA4 VUS pulls.
   - If blocked, log reason + required approvals.

4. **SpliceAI + AlphaMissense score tables**
   - Identify public BigQuery/TSV dumps for each predictor.
   - Note coordinate systems, join keys, and any licensing constraints.

5. **Protein domain / structural metadata**
   - Source: UniProt, ClinGen, published ABCA4 domain papers.
   - Record URLs/DOIs and extraction plan (manual table vs script).

6. **DNA foundation model weights**
   - List download commands for Enformer (HuggingFace), HyenaDNA, or other models we will probe.
   - Confirm size, format (PyTorch, safetensors), and GPU requirements.

7. **Storage targets**
   - Define where each asset lands locally (`data_raw/`, `external/`), expected file naming, and version tagging.

For every item above, the agent should fill a table with: `dataset | source URL | auth needed? | download cmd | size/checksum | last verified date`.

---

## 3. Package & Tooling Survey

Ask the agent to run the following checks (recording date + command output snippet):

1. **Core genomics IO**: `cyvcf2`, `pysam`, `pyensembl`, `biopython`.
2. **Annotation / modeling**: `vep` (CLI), `alphamissense` data loader, `spliceai` CLI or python pkg, `enformer-pytorch`, `evo2`, `hyenadna`.
3. **Data wrangling & infra**: `pandas`, `polars`, `numpy`, `pyarrow`, `duckdb`.
4. **Optimization / Strand SDK**: `strand-sdk` (local path + version), `hydra-core`, `mlflow`, `rich`, `pydantic`, `torch`.
5. **Workflow/orchestration**: evaluate `uv` or `rye` for dependency management plus `prefect` or `dagster` for pipelines; capture pros/cons vs sticking with `invoke`/`make`.

Verification steps:
- Use `pip index versions <pkg>` or `pip install <pkg>==` failure to find latest published version.
- Note compatibility constraints (e.g., `pyarrow>=16` due to pandas 2.2, CUDA requirements for torch 2.5+).
- Produce a recommendation table: `package | latest stable | reason to adopt | pin suggestion`.

---

## 4. Orchestration Guidance

Give the agent these prompts to answer in their report:

1. **Environment management**: Should we switch to `uv` (fast lockfile) or keep `pip-tools`? Provide steps for creating an isolated env (`uv venv`, `uv pip sync`, or `poetry install`).
2. **Data → feature → optimization flow**: propose either a lightweight `make` + `pyproject` scripts or a `prefect` flow. Include which scripts run sequentially (download → annotate → feature compute → Strand run → report render).
3. **Caching & provenance**: recommend how to store model inputs/outputs (e.g., `data_processed/parquet` + MLflow artifacts). Specify logging hooks we should wire into `strand/logging/mlflow_tracker.py`.
4. **Automation hooks**: outline CI-style checks (lint, datatype validation) and how to integrate them into `EXECUTION_PLAN.md` milestones.

---

## 5. Handoff Instructions

- Place findings in `docs/research/abca4_campaign/data_source_and_orchestration_report.md`.
- Update `EXECUTION_PLAN.md` “Data & Infra” section with a one-line summary + link.
- Flag any blockers (licensing, missing APIs) as GitHub issues or in a `BLOCKERS.md` file under this folder.

This plan should give the next agent everything needed to verify real-world data sources, capture the freshest package info, and describe how those tools orchestrate the ABCA4 campaign end-to-end.
