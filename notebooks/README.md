# ABCA4 Interactive Marimo Notebooks

This directory contains three comprehensive, interconnected Marimo notebooks that implement the complete ABCA4 variant intelligence pipeline—from raw data exploration through optimization and reporting.

## Notebooks Overview

### 1. `01_data_exploration.py` – Data Ingest & Annotation

**Purpose:** Load, filter, and annotate raw variants with functional predictions and conservation scores.

**Steps Covered:**
- **Step 0:** Scope & framing (ABCA4, Stargardt disease context, campaign goals)
- **Step 1:** Data ingest (ClinVar/gnomAD downloads, optional TSV upload, live filtering)
- **Step 2:** Annotation & deterministic features (VEP, gnomAD join, conservation, domain mapping)

**Key Features:**
- Interactive parameter cell: gene symbol, transcript, panel size, budget narrative
- File upload widget for partner TSV data
- Live filter panels: clinical significance, allele frequency, domain
- Completeness metrics for each annotation type
- Auto-export of `abca4_vus_annotated.parquet`

**Output:** `data_processed/annotations/abca4_vus_annotated.parquet`

---

### 2. `02_feature_engineering.py` – Scoring & Clustering

**Purpose:** Compute model scores (AlphaMissense, SpliceAI, LoF priors) and construct impact metrics.

**Steps Covered:**
- **Step 3:** Main model scoring (load annotated variants, add scores, visualize distributions)
- **Step 4:** Impact score construction (hand-mix weights vs. logistic regression, calibration plots)
- **Step 5:** Clustering & coverage targets (domain-based or consequence-based clusters, tau_j thresholds)

**Key Features:**
- Dual scoring modes: manual weight tuning (sliders) or logistic regression training
- Real-time calibration plots: pathogenic vs. benign score distributions
- Interactive clustering strategy selection
- Per-cluster coverage target calculation
- Auto-export of `variants_scored.parquet` with cluster assignments

**Output:** `data_processed/features/variants_scored.parquet`

---

### 3. `03_optimization_dashboard.py` – Strand Search & Reporting

**Purpose:** Run Strand optimization to select top K variants, map to experimental assays, and generate reports.

**Steps Covered:**
- **Step 6:** Strand environment & search (widget controls for K, iterations, strategy, reward weights)
- **Step 7:** Experimental mapping (consequence-to-mechanism, suggested assays, editable rationale)
- **Step 8:** Report preview (assembled Markdown report with context, approach, selected variants, assay plan)

**Key Features:**
- Interactive reward weight sliders + strategy selector (CEM/GA/Random)
- Run button for in-place Strand optimization with MLflow logging
- Consequence-to-mechanism and consequence-to-assay lookup tables
- Editable variant rationale fields
- Auto-generated Markdown report with all metadata
- Exports: CSV, JSON, and Markdown report snapshots

**Outputs:**
- `data_processed/reports/variants_selected.csv`
- `data_processed/reports/variants_selected.json`
- `data_processed/reports/report_snapshot.md`

---

### 4. `04_fasta_exploration.py` – FASTA Motif Viewer (Optional)

**Purpose:** Explore ABCA4 protein sequence motifs for ATP-binding sites and N-glycosylation patterns.

**Required Input:** `data_raw/sequences/ABCA4_P78363.fasta` (UniProt FASTA format)

**Features:**
- ATP-binding P-loop/Walker A motifs: `GxxxxGK[ST]`
- N-glycosylation motifs: `N[^P][ST][^P]`
- Interactive motif position display

**Note:** This is a utility notebook for protein sequence analysis. The required FASTA file is not included in the repository and must be obtained separately from UniProt (P78363) if motif analysis is needed.

---

## Usage

### Interactive Editing (Recommended for Development)

```bash
# Run all notebooks in interactive edit mode
marimo edit notebooks/01_data_exploration.py
marimo edit notebooks/02_feature_engineering.py
marimo edit notebooks/03_optimization_dashboard.py
```

In edit mode, cells automatically re-run when you modify them, and all widgets are interactive.

### Dashboard Mode (for Sharing / Presentation)

```bash
# Deploy as standalone dashboards
marimo run notebooks/01_data_exploration.py
marimo run notebooks/02_feature_engineering.py
marimo run notebooks/03_optimization_dashboard.py
```

Each will start an interactive web app on localhost:3000 (with port increments for multiple instances).

### Script Mode (for Programmatic Execution)

```bash
# Run as Python scripts with default parameters
python notebooks/01_data_exploration.py
python notebooks/02_feature_engineering.py
python notebooks/03_optimization_dashboard.py
```

## Data Flow

```
data_raw/
  ├── clinvar/
  ├── gnomad/
  ├── spliceai/
  └── alphamissense/
        ↓
01_data_exploration.py
        ↓
data_processed/annotations/
  └── abca4_vus_annotated.parquet
        ↓
02_feature_engineering.py
        ↓
data_processed/features/
  ├── variants_features_raw.parquet
  └── variants_scored.parquet
        ↓
03_optimization_dashboard.py
        ↓
data_processed/reports/
  ├── variants_selected.csv
  ├── variants_selected.json
  └── report_snapshot.md
```

## Configuration

### Marimo Runtime Settings (.marimo.toml)

The campaign's `.marimo.toml` configures:
- **Theme:** Light (optimized for data viz)
- **Runtime:** Lazy evaluation (cells run only when needed)
- **Package Manager:** uv (fast dependency resolution)
- **Formatting:** Auto-format with Ruff on save

### Parameter Cells

Each notebook has a dedicated **parameter cell** near the top:

**01_data_exploration.py:**
```python
gene_symbol = "ABCA4"
transcript_id = "ENST00000370225"
k_variants = 30
budget_narrative = "..."
```

**02_feature_engineering.py:**
```python
scoring_mode = "hand-mix"  # or "logistic"
clustering_mode = "domain"  # or "consequence", "manual"
```

**03_optimization_dashboard.py:**
```python
strategy = "CEM"  # or "GA", "Random"
k_variants = 30
num_iterations = 1000
```

All downstream cells automatically depend on these parameters—modify them once, and the entire pipeline updates.

## Reactivity & Data Binding

Marimo's reactive execution model means:

1. **Cell dependencies are explicit:** Each cell declares its inputs (via function parameters).
2. **Updates propagate automatically:** Change a slider, and all dependent cells re-run.
3. **No kernel state issues:** Every run starts fresh; reproducible every time.
4. **Interactive widgets:** Sliders, buttons, text inputs, dropdowns, and file uploads are first-class.

Example:
```python
@app.cell
def __(mo):
    K = mo.ui.slider(10, 200, value=30, label="Panel Size")
    return K

@app.cell
def __(K):
    # This cell re-runs whenever K changes
    selected = df[df["rank"] <= K.value]
    return selected
```

## Dependencies

Install campaign-specific requirements:

```bash
pip install -r requirements.txt
```

Key dependencies:
- **marimo** (≥0.17.8): Interactive notebooks
- **pandas, numpy, scipy:** Data manipulation & analysis
- **pyensembl:** Transcript annotation
- **pysam:** VCF parsing
- **requests:** API calls (VEP, gnomAD)
- **plotly:** Interactive visualizations

## Advanced Usage

### Connecting to External Data

Edit the `DATA_RAW_DIR` paths in each notebook to point to your data:

```python
DATA_RAW_DIR = Path("data_raw")  # default
# or
DATA_RAW_DIR = Path("/mnt/shared/gnomad_data")  # custom path
```

### Custom Scoring Functions

The notebook implements **deterministic scoring** using hand-mix weighting or logistic regression:

**Hand-mix mode:** Combines normalized AlphaMissense, SpliceAI, phyloP conservation, and LoF prior scores using adjustable weights.

**Logistic regression mode:** Trains a classifier on labeled pathogenic/benign variants to learn optimal feature combinations.

**DEMO_MODE_ENABLED:** Set `DEMO_MODE_ENABLED = True` in `src/config.py` to bypass errors with synthetic data for demonstration purposes.

For production use, ensure all required features are available. The system will fail fast with clear error messages if scoring cannot proceed.

### Custom Reward Blocks

In `03_optimization_dashboard.py`, extend the reward weight section:

```python
mo.md("### Additional Reward Blocks")
tfbs_weight = mo.ui.slider(0, 1, value=0.1, label="TFBS Conservation")
# Then add to the normalized_weights dict
```

## Troubleshooting

**Issue:** Notebooks won't load / "Missing data_processed/"

**Solution:** Create the directory structure:
```bash
mkdir -p data_processed/{annotations,features,reports}
```

**Issue:** "Module not found: pyensembl"

**Solution:** Install dependencies:
```bash
pip install -r requirements.txt
```

**Issue:** Widgets don't update downstream cells

**Solution:** Ensure you're calling `.value` on widget objects:
```python
k_value = k_slider.value  # ✓ Correct
# NOT:
k_value = k_slider  # ✗ Will cause issues
```

## Extending the Pipeline

To add a new analysis step:

1. **Create a new cell** with a meaningful name:
   ```python
   @app.cell
   def __(mo, df_input):
       """My new analysis step."""
       # Your code here
       return result
   ```

2. **Declare dependencies** via function parameters:
   ```python
   def __(mo, df_input, my_param):  # Marimo sees these dependencies
   ```

3. **Return outputs** to make them available downstream:
   ```python
   return new_dataframe, summary_dict
   ```

4. **Test reactivity:** Modify a widget that feeds into your cell, and verify it updates.

## References

- [Marimo Docs](https://docs.marimo.io/) – Complete framework reference
- [Marimo GitHub](https://github.com/marimo-team/marimo) – Source & examples
- [ABCA4 Research](./../../docs/research/) – Campaign-specific literature
- [Strand SDK](./../../README.md) – Main framework documentation

## Next Steps

Once you've run all three notebooks:

1. Review the generated `data_processed/reports/report_snapshot.md`
2. Export selected variants and share with experimental colleagues
3. Run formal optimization via `src/reward/run_abca4_optimization.py` for batch processing
4. Iterate: adjust weights, try new clustering strategies, and re-run from notebook widgets

---

