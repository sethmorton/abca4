#!/usr/bin/env python3
"""
ABCA4 Campaign - Invoke Tasks

Orchestration tasks for the ABCA4 variant triage campaign.
Run with: invoke <task-name>
"""

import os
from pathlib import Path
from invoke import task, Collection
import sys

CAMPAIGN_ROOT = Path(__file__).resolve().parent
REPO_ROOT = CAMPAIGN_ROOT.parents[1]
sys.path.insert(0, str(REPO_ROOT))

DATA_RAW = CAMPAIGN_ROOT / "data_raw"
DATA_PROCESSED = CAMPAIGN_ROOT / "data_processed"
NOTEBOOKS = CAMPAIGN_ROOT / "notebooks"
SRC = CAMPAIGN_ROOT / "src"


def _run_script(c, script_path: Path, description: str, optional: bool = False):
    if not script_path.exists():
        if optional:
            print(f"⚠️  Skipping {description} (missing {script_path})")
            return
        raise FileNotFoundError(f"Required script not found: {script_path}")
    c.run(f"uv run python {script_path}")

@task
def setup_dev(c):
    """Set up development environment with uv"""
    print("🚀 Setting up ABCA4 campaign development environment...")

    # Install uv if not present
    if not c.run("which uv", warn=True, hide=True).ok:
        print("📦 Installing uv...")
        c.run("curl -LsSf https://astral.sh/uv/install.sh | sh")

    with c.cd(str(REPO_ROOT)):
        print("🏗️ Creating virtual environment...")
        c.run("uv venv")

        print("📚 Installing dependencies...")
        c.run("uv pip sync requirements-dev.txt pyproject.toml")

        print("🔧 Installing project in development mode...")
        c.run("uv pip install -e .")

        print("🎨 Installing interactive dependencies (Marimo)...")
        c.run("uv pip install -e .[interactive]")

        extra_reqs = CAMPAIGN_ROOT / "requirements.txt"
        if extra_reqs.exists():
            print("➕ Installing campaign-specific requirements...")
            c.run(f"uv pip install -r {extra_reqs}")

    print("✅ Development environment ready!")
    print("\n🎯 Quick start:")
    print("  uv run marimo edit notebooks/01_data_exploration.py")
    print("  invoke download-data")
    print("  invoke run-pipeline")

@task
def download_data(c):
    """Download all required datasets"""
    print("📥 Downloading ABCA4 campaign datasets...")

    # Ensure data directories exist
    DATA_RAW.mkdir(exist_ok=True)
    (DATA_RAW / "clinvar").mkdir(exist_ok=True)
    (DATA_RAW / "gnomad").mkdir(exist_ok=True)
    (DATA_RAW / "spliceai").mkdir(exist_ok=True)
    (DATA_RAW / "alphamissense").mkdir(exist_ok=True)

    with c.cd(str(REPO_ROOT)):
        print("🧬 Downloading ClinVar data...")
        _run_script(c, SRC / "data" / "download_clinvar.py", "ClinVar download")

        print("🧬 Downloading gnomAD data...")
        _run_script(c, SRC / "data" / "download_gnomad.py", "gnomAD download")

        print("🧬 Downloading SpliceAI data...")
        _run_script(c, SRC / "data" / "download_spliceai.py", "SpliceAI download")

        print("🧬 Downloading AlphaMissense data...")
        _run_script(c, SRC / "data" / "download_alphamissense.py", "AlphaMissense download")

    print("✅ All datasets downloaded!")

@task
def process_variants(c):
    """Process and filter ABCA4 variants from raw data"""
    print("🔍 Processing ABCA4 variants...")

    DATA_PROCESSED.mkdir(exist_ok=True)
    (DATA_PROCESSED / "variants").mkdir(exist_ok=True)

    with c.cd(str(REPO_ROOT)):
        _run_script(c, SRC / "data" / "filter_abca4_variants.py", "variant filtering")

    print("✅ ABCA4 variants processed!")

@task
def annotate_variants(c):
    """Add transcript and functional annotations to variants"""
    print("📝 Annotating variants...")

    (DATA_PROCESSED / "annotations").mkdir(exist_ok=True)

    with c.cd(str(REPO_ROOT)):
        _run_script(c, SRC / "annotation" / "annotate_transcripts.py", "variant annotation")

    print("✅ Variants annotated!")

@task
def compute_features(c):
    """Compute all feature matrices for variants"""
    print("🧮 Computing features...")

    (DATA_PROCESSED / "features").mkdir(exist_ok=True)

    with c.cd(str(REPO_ROOT)):
        _run_script(c, SRC / "features" / "missense.py", "missense feature computation")
        _run_script(c, SRC / "features" / "splice.py", "splice feature computation")
        _run_script(c, SRC / "features" / "regulatory.py", "regulatory features")
        _run_script(c, SRC / "features" / "conservation.py", "conservation features")
        _run_script(c, SRC / "features" / "assemble_features.py", "feature assembly")

    print("✅ Features computed!")

@task
def run_optimization(c):
    """Run Strand optimization campaign"""
    print("🚀 Running optimization campaign...")
    with c.cd(str(REPO_ROOT)):
        _run_script(
            c,
            SRC / "reward" / "run_abca4_optimization.py",
            "optimization runner",
        )

    print("✅ Optimization complete!")

@task
def generate_report(c):
    """Generate final reports and dashboards"""
    print("📊 Generating reports...")

    (DATA_PROCESSED / "reports").mkdir(exist_ok=True)

    with c.cd(str(REPO_ROOT)):
        _run_script(c, SRC / "reporting" / "generate_snapshot.py", "report generation")

    print("✅ Reports generated!")

@task
def run_pipeline(c):
    """Run the complete ABCA4 pipeline end-to-end"""
    print("🔬 Running complete ABCA4 campaign pipeline...")

    download_data(c)
    process_variants(c)
    annotate_variants(c)
    compute_features(c)
    run_optimization(c)
    generate_report(c)

    print("🎉 ABCA4 campaign complete!")

@task
def explore_data(c):
    """Launch interactive data exploration notebook"""
    print("🔬 Launching data exploration notebook...")
    with c.cd(str(REPO_ROOT)):
        if (NOTEBOOKS / "01_data_exploration.py").exists():
            c.run(f"uv run marimo edit {NOTEBOOKS}/01_data_exploration.py")
        else:
            print("⚠️  Data exploration notebook missing; skipping")

@task
def tune_features(c):
    """Launch interactive feature engineering notebook"""
    print("🔧 Launching feature engineering notebook...")
    with c.cd(str(REPO_ROOT)):
        if (NOTEBOOKS / "02_feature_engineering.py").exists():
            c.run(f"uv run marimo edit {NOTEBOOKS}/02_feature_engineering.py")
        else:
            print("⚠️  Feature engineering notebook missing; skipping")

@task
def optimize_interactive(c):
    """Launch interactive optimization dashboard"""
    print("🚀 Launching optimization dashboard...")
    with c.cd(str(REPO_ROOT)):
        if (NOTEBOOKS / "03_optimization_dashboard.py").exists():
            c.run(f"uv run marimo run {NOTEBOOKS}/03_optimization_dashboard.py")
        else:
            print("⚠️  Optimization notebook missing; skipping")

@task
def clean_data(c):
    """Clean intermediate data files"""
    print("🧹 Cleaning intermediate data...")

    import shutil

    # Remove processed data but keep raw downloads
    if DATA_PROCESSED.exists():
        shutil.rmtree(DATA_PROCESSED)
        print("✅ Processed data cleaned!")

@task
def test_pipeline(c):
    """Run tests for the ABCA4 pipeline"""
    print("🧪 Running pipeline tests...")
    tests_path = REPO_ROOT / "tests" / "abca4"
    if not tests_path.exists():
        print("⚠️  No ABCA4-specific tests found; skipping")
        return

    with c.cd(str(REPO_ROOT)):
        c.run(f"uv run pytest {tests_path} -v")

    print("✅ Tests complete!")

@task
def cro_parse(c):
    """Parse variant report into CRO variant panel (Stage 1)"""
    print("🔍 Parsing variant report into CRO variant panel...")
    c.run(f"cd {CAMPAIGN_ROOT} && uv run python src/cro/parser.py")
    print("✅ CRO variant panel parsed!")

@task
def cro_annotate(c):
    """Annotate variants with mechanisms (Stage 2)"""
    print("🧬 Annotating variants with molecular mechanisms...")
    c.run(f"cd {CAMPAIGN_ROOT} && uv run python src/cro/mechanism.py")
    print("✅ Variants annotated with mechanisms!")

@task
def cro_assign(c):
    """Assign assay modules to variants (Stage 3)"""
    print("🧪 Assigning assay modules to variants...")
    c.run(f"cd {CAMPAIGN_ROOT} && uv run python src/cro/assay_mapper.py")
    print("✅ Assay modules assigned!")

@task
def cro_workpackages(c):
    """Generate work packages (Stage 4)"""
    print("📦 Generating work packages...")
    c.run(f"cd {CAMPAIGN_ROOT} && uv run python src/cro/workpackages.py")
    print("✅ Work packages generated!")

@task
def cro_designs(c):
    """Generate experimental designs (Stage 5)"""
    print("📊 Generating experimental designs...")
    c.run(f"cd {CAMPAIGN_ROOT} && uv run python src/cro/designs.py")
    print("✅ Experimental designs generated!")

@task
def cro_deliverables(c):
    """Generate deliverable specifications (Stage 6)"""
    print("📋 Generating deliverable specifications...")
    c.run(f"cd {CAMPAIGN_ROOT} && uv run python src/cro/deliverables.py")
    print("✅ Deliverable specifications generated!")

@task
def cro_validate(c):
    """Validate CRO pipeline outputs"""
    print("🔍 Validating CRO pipeline...")
    c.run(f"cd {CAMPAIGN_ROOT} && uv run python src/cro/cro_validate.py")
    print("✅ CRO validation complete!")

@task
def cro_plan(c):
    """Generate complete CRO study plan (Stages 1-7)"""
    print("📄 Generating complete CRO study plan...")

    # Run all stages in sequence
    cro_parse(c)
    cro_annotate(c)
    cro_assign(c)
    cro_workpackages(c)
    cro_designs(c)
    cro_deliverables(c)

    # Run validation
    print("🔍 Running final validation...")
    cro_validate(c)

    # Generate final study plan
    print("🎨 Rendering final CRO study plan...")
    c.run(f"cd {CAMPAIGN_ROOT} && uv run python src/reporting/generate_cro_plan.py")

    print("🎉 CRO study plan complete!")
    print("\n📂 Outputs available in:")
    print("  data_processed/cro/ - All intermediate artifacts")
    print("  data_processed/cro/validation_report.json - Validation results")
    print("  data_processed/reports/cro_study_plan.md - Final study plan")

@task
def cro_dashboard(c):
    """Launch CRO planning dashboard"""
    print("🚀 Launching CRO planning dashboard...")
    if (NOTEBOOKS / "05_cro_plan.py").exists():
        c.run(f"cd {CAMPAIGN_ROOT} && uv run marimo run {NOTEBOOKS}/05_cro_plan.py")
    else:
        print("⚠️  CRO dashboard notebook missing; run 'invoke cro.plan' first")

# Create collections for organization
data_ns = Collection('data')
data_ns.add_task(download_data, 'download')
data_ns.add_task(process_variants, 'process')

features_ns = Collection('features')
features_ns.add_task(compute_features, 'compute')

cro_ns = Collection('cro')
cro_ns.add_task(cro_parse, 'parse')
cro_ns.add_task(cro_annotate, 'annotate')
cro_ns.add_task(cro_assign, 'assign')
cro_ns.add_task(cro_workpackages, 'workpackages')
cro_ns.add_task(cro_designs, 'designs')
cro_ns.add_task(cro_deliverables, 'deliverables')
cro_ns.add_task(cro_validate, 'validate')
cro_ns.add_task(cro_plan, 'plan')
cro_ns.add_task(cro_dashboard, 'dashboard')

notebook_ns = Collection('notebook')
notebook_ns.add_task(explore_data, 'explore')
notebook_ns.add_task(tune_features, 'tune')
notebook_ns.add_task(optimize_interactive, 'optimize')

# Main namespace
ns = Collection()
ns.add_task(setup_dev)
ns.add_task(run_pipeline)
ns.add_task(run_optimization)
ns.add_task(generate_report)
ns.add_task(clean_data)
ns.add_task(test_pipeline)
ns.add_collection(data_ns)
ns.add_collection(features_ns)
ns.add_collection(cro_ns)
ns.add_collection(notebook_ns)
