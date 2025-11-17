"""Generate CRO Study Plan Document

Renders comprehensive study plan from all CRO pipeline stages using Jinja templates.
"""

from __future__ import annotations

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict

import pandas as pd
from jinja2 import Environment, FileSystemLoader

try:
    # When run as module
    from ..cro.assay_mapper import load_assay_catalog
    from ..cro.cro_types import (
        VariantPanelRowData, WorkPackageData, DesignSummaryData, DeliverableSpecData
    )
    from ..cro.deliverables import DeliverableSpecs
    from ..cro.designs import DesignSummaries
    from ..cro.parser import VariantPanel
    from ..cro.workpackages import WorkPackages
except ImportError:
    # When run as standalone script
    import sys
    from pathlib import Path
    # Add the src directory to sys.path
    src_dir = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(src_dir))
    from cro.assay_mapper import load_assay_catalog
    from cro.cro_types import (
        VariantPanelRowData, WorkPackageData, DesignSummaryData, DeliverableSpecData
    )
    from cro.deliverables import DeliverableSpecs
    from cro.designs import DesignSummaries
    from cro.parser import VariantPanel
    from cro.workpackages import WorkPackages

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(name)s | %(levelname)s | %(message)s",
)
LOGGER = logging.getLogger(__name__)

CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]
CRO_DATA_DIR = CAMPAIGN_ROOT / "data_processed" / "cro"
REPORTS_DIR = CAMPAIGN_ROOT / "data_processed" / "reports"
TEMPLATES_DIR = CAMPAIGN_ROOT / "src" / "reporting" / "templates"


def load_all_cro_data() -> tuple[list[VariantPanelRowData], list[WorkPackageData], list[DesignSummaryData], list[DeliverableSpecData]]:
    """Load all CRO pipeline data with type validation."""
    # Load variant panel
    variant_df = pd.read_parquet(CRO_DATA_DIR / "variant_panel.parquet")
    variant_panel: list[VariantPanelRowData] = []
    for _, row in variant_df.iterrows():
        # Convert to typed dict - pandas dtypes need to be converted to Python types
        variant_data = {
            'variant_id': str(row['variant_id']),
            'gene': str(row['gene']),
            'rank': int(row['rank']),
            'chrom': str(row['chrom']),
            'pos': int(row['pos']),
            'ref': str(row['ref']),
            'alt': str(row['alt']),
            'consequence': str(row['consequence']),
            'domain': str(row['domain']),
            'model_score': float(row['model_score']),
            'impact_score': float(row['impact_score']),
            'cluster_id': str(row['cluster_id']),
            'coverage': float(row['coverage']),
            'assay_hint': str(row['assay_hint']),
            'selection_strategy': str(row['selection_strategy']),
            'panel_size': int(row['panel_size']),
            'panel_id': str(row['panel_id']),
            'source_md_path': str(row['source_md_path'])
        }
        variant_panel.append(variant_data)

    # Load work packages with validation
    work_packages: list[WorkPackageData] = []
    with open(CRO_DATA_DIR / "work_packages.jsonl", "r") as f:
        for line in f:
            if line.strip():
                wp_data = json.loads(line)
                # Basic validation - ensure required fields exist
                required_wp_fields = ['wp_id', 'variant_ids', 'assay_module', 'gene', 'objective']
                for field in required_wp_fields:
                    if field not in wp_data:
                        raise ValueError(f"Work package missing required field: {field}")
                work_packages.append(wp_data)

    # Load design summaries with validation
    with open(CRO_DATA_DIR / "design_summaries.json", "r") as f:
        design_summaries_raw = json.load(f)
    design_summaries: list[DesignSummaryData] = []
    for ds_data in design_summaries_raw:
        # Basic validation
        required_ds_fields = ['wp_id', 'factors', 'design_type', 'tech_reps', 'bio_reps', 'controls']
        for field in required_ds_fields:
            if field not in ds_data:
                raise ValueError(f"Design summary missing required field: {field}")
        design_summaries.append(ds_data)

    # Load deliverable specs with validation
    with open(CRO_DATA_DIR / "deliverable_specs.json", "r") as f:
        deliverable_specs_raw = json.load(f)
    deliverable_specs: list[DeliverableSpecData] = []
    for ds_data in deliverable_specs_raw:
        # Basic validation
        required_ds_fields = ['wp_id', 'primary_metrics', 'raw_returns', 'summary_columns', 'qc_expectations']
        for field in required_ds_fields:
            if field not in ds_data:
                raise ValueError(f"Deliverable spec missing required field: {field}")
        deliverable_specs.append(ds_data)

    return variant_panel, work_packages, design_summaries, deliverable_specs


def _create_variant_to_workpackage_mapping(
    variant_panel: list[VariantPanelRowData],
    work_packages: list[WorkPackageData]
) -> Dict[str, list]:
    """Create mapping from variants to work packages for appendix."""
    mapping = {}

    # Create lookup from wp_id to wp
    wp_lookup = {wp["wp_id"]: wp for wp in work_packages}

    for wp in work_packages:
        for variant_id in wp["variant_ids"]:
            if variant_id not in mapping:
                mapping[variant_id] = []
            mapping[variant_id].append({
                "wp_id": wp["wp_id"],
                "assay_module": wp["assay_module"],
                "gene": wp["gene"],
            })

    return mapping


def _calculate_study_summary_stats(
    variant_panel: list[VariantPanelRowData],
    work_packages: list[WorkPackageData],
    design_summaries: list[DesignSummaryData]
) -> Dict[str, str | int | float | list]:
    """Calculate summary statistics for the study plan."""
    stats = {
        "total_variants": len(variant_panel),
        "total_work_packages": len(work_packages),
        "assay_modules_used": list(set(wp["assay_module"] for wp in work_packages)),
        "genes_studied": list(set(vp["gene"] for vp in variant_panel)),
        "estimated_total_samples": 0,
        "estimated_duration_weeks": 0,
    }

    # Calculate total samples and duration
    assay_catalog = load_assay_catalog()
    for wp in work_packages:
        module = assay_catalog[wp["assay_module"]]
        variant_count = len(wp["variant_ids"])

        # Estimate samples (variants + controls) × replicates
        ds = next((d for d in design_summaries if d["wp_id"] == wp["wp_id"]), None)
        if ds:
            samples_per_condition = ds["tech_reps"] * ds["bio_reps"]
            total_conditions = variant_count + len(ds["controls"])
            stats["estimated_total_samples"] += total_conditions * samples_per_condition

        # Add duration
        stats["estimated_duration_weeks"] = max(
            stats["estimated_duration_weeks"],
            module.estimated_duration_weeks
        )

    return stats


def render_study_plan(
    variant_panel: list,
    work_packages: list,
    design_summaries: list,
    deliverable_specs: list,
    output_path: Path = REPORTS_DIR / "cro_study_plan.md"
) -> Path:
    """Render the complete CRO study plan using Jinja template.

    Args:
        variant_panel: Variant panel data
        work_packages: Work package data
        design_summaries: Design summary data
        deliverable_specs: Deliverable specification data
        output_path: Path to save the rendered markdown

    Returns:
        Path to the rendered file
    """
    LOGGER.info("Rendering CRO study plan to %s", output_path)

    # Set up Jinja environment
    templates_dir = TEMPLATES_DIR
    templates_dir.mkdir(parents=True, exist_ok=True)

    env = Environment(loader=FileSystemLoader(templates_dir))
    template = env.get_template("cro_study_plan.md.jinja")

    # Prepare template context
    context = {
        "generation_date": datetime.now().strftime("%Y-%m-%d"),
        "variant_panel": variant_panel,
        "work_packages": work_packages,
        "design_summaries": design_summaries,
        "deliverable_specs": deliverable_specs,
        "variant_to_wp_mapping": _create_variant_to_workpackage_mapping(variant_panel, work_packages),
        "study_stats": _calculate_study_summary_stats(variant_panel, work_packages, design_summaries),
    }

    # Render template
    rendered_content = template.render(**context)

    # Save to file
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(rendered_content)

    LOGGER.info("Successfully rendered CRO study plan to %s", output_path)
    return output_path


def create_template() -> Path:
    """Create the Jinja template if it doesn't exist."""
    template_path = TEMPLATES_DIR / "cro_study_plan.md.jinja"
    if template_path.exists():
        return template_path

    template_content = """# CRO Study Plan: ABCA4 Variant Functional Validation

**Generated:** {{ generation_date }}
**Campaign:** ABCA4 Variant Triage
**Total Variants:** {{ study_stats.total_variants }}
**Work Packages:** {{ study_stats.total_work_packages }}

---

## Executive Summary

This study plan outlines a comprehensive functional validation campaign for {{ study_stats.total_variants }} ABCA4 variants identified through computational prioritization. The campaign is organized into {{ study_stats.total_work_packages }} work packages spanning {{ study_stats.assay_modules_used | length }} assay modalities.

### Study Overview
- **Genes:** {{ study_stats.genes_studied | join(", ") }}
- **Assay Modules:** {{ study_stats.assay_modules_used | join(", ") }}
- **Estimated Samples:** {{ study_stats.estimated_total_samples | default("TBD") }}
- **Estimated Duration:** {{ study_stats.estimated_duration_weeks | default("TBD") }} weeks

---

## Work Packages

{% for wp in work_packages %}
### {{ wp.wp_id }}

**Objective:** {{ wp.objective }}

**Assay Module:** {{ wp.assay_module }}
**Variants:** {{ wp.variant_ids | length }}
**Gene:** {{ wp.gene }}

#### Materials Provided
{% for material in wp.materials_provided %}
- {{ material }}
{% endfor %}

#### Materials Needed
{% for material in wp.materials_needed %}
- {{ material }}
{% endfor %}

#### Experimental Design
{% set ds = design_summaries | selectattr("wp_id", "equalto", wp.wp_id) | first %}
- **Design Type:** {{ ds.design_type }}
- **Technical Replicates:** {{ ds.tech_reps }}
- **Biological Replicates:** {{ ds.bio_reps }}
- **Factors:**
{% for factor in ds.factors %}
  - {{ factor.name }}: {{ factor.levels | join(", ") }}
{% endfor %}
- **Controls:** {{ ds.controls | join(", ") }}

#### Deliverables
{% set ds_spec = deliverable_specs | selectattr("wp_id", "equalto", wp.wp_id) | first %}
**Primary Metrics:**
{% for metric in ds_spec.primary_metrics %}
- {{ metric }}
{% endfor %}

**Raw Data Returns:**
{% for raw_return in ds_spec.raw_returns %}
- {{ raw_return }}
{% endfor %}

**Summary Data Columns:**
{% for column in ds_spec.summary_columns %}
- {{ column }}
{% endfor %}

**QC Expectations:**
{% for qc in ds_spec.qc_expectations %}
- {{ qc }}
{% endfor %}

#### CRO Notes
{{ wp.cro_notes }}

---

{% endfor %}

## Data & Reporting

### Data Submission Format
All data should be submitted in the following structure:
```
cro_results/
├── {wp_id}/
│   ├── raw_data/
│   │   └── {raw_returns}
│   ├── processed_data/
│   │   └── {wp_id}_summary.csv
│   └── reports/
│       └── {wp_id}_report.pdf
```

### Quality Control Requirements
- All assays must include appropriate positive and negative controls
- Replicate concordance must meet specifications in QC expectations
- Raw data and processing pipelines must be provided
- Any deviations from protocol must be documented

### Timeline Expectations
- Kick-off meeting within 1 week of contract execution
- Regular progress updates bi-weekly
- Preliminary data review at midpoint
- Final data package within 2 weeks of study completion

---

## Appendix: Variant to Work Package Mapping

| Variant ID | Work Packages |
|------------|---------------|
{% for variant_id, wps in variant_to_wp_mapping.items() %}
| {{ variant_id }} | {% for wp in wps %}{{ wp.wp_id }} ({{ wp.assay_module }}) {% if not loop.last %}; {% endif %}{% endfor %} |
{% endfor %}

---

*Generated automatically from ABCA4 variant prioritization pipeline*
"""

    with open(template_path, "w") as f:
        f.write(template_content)

    LOGGER.info("Created template at %s", template_path)
    return template_path


def main() -> None:
    """Generate the CRO study plan document."""
    # Create template if needed
    create_template()

    # Load all data
    variant_panel, work_packages, design_summaries, deliverable_specs = load_all_cro_data()

    # Render study plan
    output_path = REPORTS_DIR / "cro_study_plan.md"
    render_study_plan(
        variant_panel, work_packages, design_summaries, deliverable_specs, output_path
    )

    LOGGER.info("CRO study plan generation complete: %s", output_path)


if __name__ == "__main__":
    main()
