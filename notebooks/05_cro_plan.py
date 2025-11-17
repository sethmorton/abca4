#!/usr/bin/env python3
"""
ABCA4 Campaign – CRO Study Plan Dashboard

Interactive dashboard for reviewing and adjusting CRO study plans.
Allows human-in-the-loop tuning of assay assignments and experimental designs.

Run interactively:  marimo edit notebooks/05_cro_plan.py
Run as dashboard:   marimo run notebooks/05_cro_plan.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import marimo

__generated_with = "0.17.8"
app = marimo.App()


@app.cell
def __():
    """Import core libraries and setup."""
    import marimo as mo
    import pandas as pd
    import numpy as np
    import json
    import logging
    from pathlib import Path
    from typing import Dict, List, Optional
    from datetime import datetime

    # Setup logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    return mo, pd, np, json, logger, Path, Dict, List, Optional, datetime


@app.cell
def __(mo):
    """Dashboard header and navigation."""
    mo.md("""
    # 🧬 ABCA4 CRO Study Plan Dashboard

    **Interactive Review & Tuning Interface**

    This dashboard allows you to review the complete CRO study plan pipeline,
    inspect assay assignments, experimental designs, and make final adjustments
    before generating the locked study plan document.

    ---
    """)

    # Navigation tabs
    tabs = mo.ui.tabs({
        "📊 Overview": "overview",
        "🔬 Assay Assignments": "assignments",
        "📦 Work Packages": "workpackages",
        "🧪 Experimental Designs": "designs",
        "📋 Deliverables": "deliverables",
        "📄 Generate Plan": "generate"
    })

    return tabs,


@app.cell
def __(mo, pd, json, Path):
    """Load CRO pipeline data."""
    CRO_DIR = Path("data_processed/cro")

    # Load all pipeline outputs
    data = {}

    try:
        # Variant panel
        if (CRO_DIR / "variant_panel.parquet").exists():
            data["variants"] = pd.read_parquet(CRO_DIR / "variant_panel.parquet")

        # Work packages
        if (CRO_DIR / "work_packages.jsonl").exists():
            _work_packages = []
            with open(CRO_DIR / "work_packages.jsonl", "r") as f:
                for line in f:
                    if line.strip():
                        _work_packages.append(json.loads(line))
            data["work_packages"] = _work_packages

        # Assay assignments
        if (CRO_DIR / "assay_assignments.json").exists():
            with open(CRO_DIR / "assay_assignments.json", "r") as f:
                data["assignments"] = json.load(f)

        # Design summaries
        if (CRO_DIR / "design_summaries.json").exists():
            with open(CRO_DIR / "design_summaries.json", "r") as f:
                data["designs"] = json.load(f)

        # Deliverable specs
        if (CRO_DIR / "deliverable_specs.json").exists():
            with open(CRO_DIR / "deliverable_specs.json", "r") as f:
                data["deliverables"] = json.load(f)

        # Mechanism annotations
        if (CRO_DIR / "mechanism_panel.json").exists():
            with open(CRO_DIR / "mechanism_panel.json", "r") as f:
                data["mechanisms"] = json.load(f)

        loaded_successfully = True

    except Exception as e:
        mo.callout(f"Error loading CRO data: {e}", kind="danger")
        loaded_successfully = False
        data = {}

    return data, loaded_successfully, CRO_DIR


@app.cell
def __(tabs, data, loaded_successfully, mo, pd):
    """Overview dashboard."""
    if tabs.value == "overview":
        if not loaded_successfully:
            mo.callout("CRO data not found. Run `invoke cro.plan` first.", kind="warning")
        else:
            mo.md("## 📊 Study Plan Overview")

            # Summary statistics
            variants_df = data.get("variants", pd.DataFrame())
            _work_packages = data.get("work_packages", [])

            if not variants_df.empty:
                mo.md(f"""
                ### Campaign Summary
                - **Total Variants:** {len(variants_df)}
                - **Work Packages:** {len(_work_packages)}
                - **Assay Modules:** {len(set([wp['assay_module'] for wp in _work_packages]))}
                - **Genes:** {', '.join(variants_df['gene'].unique())}
                """)

            # Work package summary
            if _work_packages:
                wp_summary = pd.DataFrame([
                    {
                        "Work Package": wp["wp_id"],
                        "Assay Module": wp["assay_module"],
                        "Variants": len(wp["variant_ids"]),
                        "Duration (weeks)": wp["cro_notes"].split("weeks")[0].split()[-1] if "weeks" in wp["cro_notes"] else "TBD"
                    }
                    for wp in work_packages
                ])
                mo.table(wp_summary, selection=None)

            # Variant distribution by consequence
            if not variants_df.empty:
                consequence_counts = variants_df["consequence"].value_counts()
                mo.bar_chart(consequence_counts, title="Variants by Consequence Type")


@app.cell
def __(tabs, data, loaded_successfully, mo, pd):
    """Assay assignments review."""
    if tabs.value == "assignments":
        if not loaded_successfully:
            mo.callout("CRO data not found. Run `invoke cro.plan` first.", kind="warning")
        else:
            mo.md("## 🔬 Assay Assignments Review")

            assignments = data.get("assignments", [])
            mechanisms = data.get("mechanisms", [])

            if assignments:
                # Create summary of assignments
                assignment_summary = []
                for assignment in assignments:
                    variant_id = assignment["variant_id"]
                    assay_modules = assignment["assay_modules"]

                    # Find mechanism annotation
                    mech_annotation = next(
                        (m for m in mechanisms if m["variant_id"] == variant_id),
                        {"mechanism_tags": []}
                    )

                    assignment_summary.append({
                        "Variant": variant_id,
                        "Assay Modules": ", ".join(assay_modules),
                        "Mechanisms": ", ".join(mech_annotation["mechanism_tags"]),
                        "Module Count": len(assay_modules)
                    })

                summary_df = pd.DataFrame(assignment_summary)

                # Summary statistics
                mo.md("### Assignment Statistics")
                module_counts = {}
                for assignment in assignments:
                    for module in assignment["assay_modules"]:
                        module_counts[module] = module_counts.get(module, 0) + 1

                stats_df = pd.DataFrame([
                    {"Assay Module": module, "Variants Assigned": count}
                    for module, count in module_counts.items()
                ])
                mo.table(stats_df, selection=None)

                # Detailed assignments table
                mo.md("### Detailed Assignments")
                mo.table(summary_df, selection="multi")


@app.cell
def __(tabs, data, loaded_successfully, mo):
    """Work packages review."""
    if tabs.value == "workpackages":
        if not loaded_successfully:
            mo.callout("CRO data not found. Run `invoke cro.plan` first.", kind="warning")
        else:
            mo.md("## 📦 Work Packages Review")

            _work_packages = data.get("work_packages", [])

            if _work_packages:
                for wp in _work_packages:
                    with mo.expander(f"📋 {wp['wp_id']} - {wp['assay_module']}"):
                        mo.md(f"""
                        **Objective:** {wp['objective']}

                        **Variants:** {len(wp['variant_ids'])}
                        - {', '.join(wp['variant_ids'][:5])}{'...' if len(wp['variant_ids']) > 5 else ''}

                        **Materials Provided:**
                        {chr(10).join(f'- {item}' for item in wp['materials_provided'])}

                        **Materials Needed:**
                        {chr(10).join(f'- {item}' for item in wp['materials_needed'])}

                        **CRO Notes:** {wp['cro_notes']}
                        """)


@app.cell
def __(tabs, data, loaded_successfully, mo, pd, Path, CRO_DIR):
    """Experimental designs review."""
    if tabs.value == "designs":
        if not loaded_successfully:
            mo.callout("CRO data not found. Run `invoke cro.plan` first.", kind="warning")
        else:
            mo.md("## 🧪 Experimental Designs Review")

            designs = data.get("designs", [])

            if designs:
                for design in designs:
                    _wp_id = design["wp_id"]
                    with mo.expander(f"🧫 {_wp_id} Design"):

                        mo.md(f"""
                        **Design Type:** {design['design_type']}
                        **Technical Replicates:** {design['tech_reps']}
                        **Biological Replicates:** {design['bio_reps']}
                        **Controls:** {', '.join(design['controls'])}
                        """)

                        # Show factors
                        if design["factors"]:
                            factors_df = pd.DataFrame([
                                {"Factor": f["name"], "Levels": ", ".join(map(str, f["levels"]))}
                                for f in design["factors"]
                            ])
                            mo.table(factors_df, selection=None)

                        # Show design CSV if it exists
                        design_csv = CRO_DIR / "designs" / f"{_wp_id}_design.csv"
                        if design_csv.exists():
                            design_df = pd.read_csv(design_csv)
                            mo.md(f"**Total Conditions:** {len(design_df)}")
                            mo.table(design_df.head(10), selection=None)
                            if len(design_df) > 10:
                                mo.md(f"*... and {len(design_df) - 10} more conditions*")


@app.cell
def __(tabs, data, loaded_successfully, mo):
    """Deliverables review."""
    if tabs.value == "deliverables":
        if not loaded_successfully:
            mo.callout("CRO data not found. Run `invoke cro.plan` first.", kind="warning")
        else:
            mo.md("## 📋 Deliverables & QC Specifications")

            deliverables = data.get("deliverables", [])

            if deliverables:
                for spec in deliverables:
                    _wp_id = spec["wp_id"]
                    with mo.expander(f"📄 {_wp_id} Deliverables"):

                        mo.md(f"""
                        **Primary Metrics:**
                        {chr(10).join(f'- {metric}' for metric in spec['primary_metrics'])}

                        **Raw Data Returns:**
                        {chr(10).join(f'- {ret}' for ret in spec['raw_returns'])}

                        **Summary Data Columns:**
                        {chr(10).join(f'- {col}' for col in spec['summary_columns'])}

                        **QC Expectations:**
                        {chr(10).join(f'- {qc}' for qc in spec['qc_expectations'])}
                        """)


@app.cell
def __(tabs, mo):
    """Generate final study plan."""
    if tabs.value == "generate":
        mo.md("## 📄 Generate Final Study Plan")

        mo.md("""
        Review all sections above and make any final adjustments needed.
        When ready, generate the locked study plan document.
        """)

        mo.md("""
        ### Generate Study Plan

        Click the button below to generate the final CRO study plan document.
        This will create a comprehensive markdown document with all work packages,
        experimental designs, and CRO instructions.
        """)

        # Generate button (simplified for Marimo)
        mo.md("""
        **Command to run:**
        ```bash
        uv run invoke cro.plan
        # or
        uv run python src/reporting/generate_cro_plan.py
        ```

        **Output:** `data_processed/reports/cro_study_plan.md`
        """)

        mo.callout("ℹ️ The study plan document will be locked and ready for CRO submission.", kind="info")


@app.cell
def __(mo):
    """Footer."""
    mo.md("""
    ---
    *ABCA4 CRO Study Plan Dashboard | Generated with Marimo*
    """)


if __name__ == "__main__":
    app.run()
