"""Stage 5: Experimental Design Generator

Creates experimental designs for each work package with factors, replicates, and controls.
"""

from __future__ import annotations

import itertools
import json
import logging
import random
from pathlib import Path
from typing import Dict, List

import pandas as pd

try:
    # When run as module
    from .cro_types import (
        AssayModuleId, DesignSummaries, DesignSummary, DesignType, Factor, WorkPackage, WorkPackages
    )
    from .workpackages import save_work_packages
except ImportError:
    # When run as standalone script
    import sys
    from pathlib import Path
    # Add the src directory to sys.path
    src_dir = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(src_dir))
    from cro.cro_types import (
        AssayModuleId, DesignSummaries, DesignSummary, DesignType, Factor, WorkPackage, WorkPackages
    )
    from cro.workpackages import save_work_packages

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(name)s | %(levelname)s | %(message)s",
)
LOGGER = logging.getLogger(__name__)

CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]
CRO_DATA_DIR = CAMPAIGN_ROOT / "data_processed" / "cro"


def load_work_packages() -> WorkPackages:
    """Load work packages from Stage 4 output."""
    jsonl_path = CRO_DATA_DIR / "work_packages.jsonl"
    if not jsonl_path.exists():
        raise FileNotFoundError(f"Work packages not found: {jsonl_path}. Run Stage 4 first.")

    work_packages = []
    with open(jsonl_path, "r") as f:
        for line in f:
            if line.strip():
                data = json.loads(line)
                wp = WorkPackage(
                    wp_id=data["wp_id"],
                    gene=data["gene"],
                    assay_module=data["assay_module"],
                    objective=data["objective"],
                    variant_ids=data["variant_ids"],
                    materials_provided=data["materials_provided"],
                    materials_needed=data["materials_needed"],
                    cro_notes=data["cro_notes"],
                )
                work_packages.append(wp)

    return work_packages


def _get_design_factors(
    work_package: WorkPackage
) -> List[Factor]:
    """Generate experimental factors for a work package."""
    factors = []

    # Primary factor: variant
    variant_levels = ["WT"] + work_package.variant_ids  # WT control + variants
    factors.append(Factor(
        name="variant",
        levels=variant_levels
    ))

    # Assay-specific factors
    assay_module = work_package.assay_module

    if assay_module == "DSF_SEC":
        # Temperature gradient for DSF
        factors.append(Factor(
            name="temperature",
            levels=[25, 30, 37, 42, 47, 52, 57, 62, 67, 72]  # °C
        ))
    elif assay_module == "FUNCTIONAL":
        # Generic functional assay - substrate concentration
        factors.append(Factor(
            name="substrate_concentration",
            levels=[0.1, 1.0, 10.0]  # µM
        ))
    elif assay_module == "TRAFFICKING":
        # Time course for trafficking
        factors.append(Factor(
            name="time_hours",
            levels=[4, 8, 12, 24, 48]  # hours post-transfection
        ))
    elif assay_module == "SPLICE_MINIGENE":
        # No additional factors - binary splicing outcome
        pass
    elif assay_module == "RNA_SEQ_REPORTER":
        # Treatment conditions
        factors.append(Factor(
            name="treatment",
            levels=["untreated", "treated"]
        ))
    elif assay_module == "TRANSCRIPTIONAL_REPORTER":
        # Promoter induction levels
        factors.append(Factor(
            name="inducer_concentration",
            levels=[0, 0.1, 1.0, 10.0]  # µM
        ))

    return factors


def _get_design_type(factors: List[Factor], variant_count: int) -> DesignType:
    """Determine the experimental design type with strict validation."""
    factor_count = len(factors)

    if factor_count == 1:
        design_type = "one_factor"
    elif factor_count == 2 and variant_count <= 10:
        design_type = "full_factorial"
    else:
        design_type = "fractional"  # For complex designs with many variants

    # Validate the design type is allowed
    valid_types = ["full_factorial", "one_factor", "fractional"]
    if design_type not in valid_types:
        raise ValueError(
            f"❌ INVALID DESIGN TYPE: {design_type}\n"
            f"Allowed types: {valid_types}\n"
            f"Factors: {factor_count}, Variants: {variant_count}"
        )

    return design_type


def _get_replicates(work_package: WorkPackage) -> tuple[int, int]:
    """Determine technical and biological replicate counts."""
    variant_count = len(work_package.variant_ids)
    assay_module = work_package.assay_module

    # Base replicates
    tech_reps = 3  # Technical replicates
    bio_reps = 2   # Biological replicates

    # Adjust based on assay type and variant count
    if assay_module in ["RNA_SEQ_REPORTER"]:
        # High-precision assays need more replicates
        tech_reps = 4
        bio_reps = 3
    elif assay_module == "DSF_SEC":
        # Biophysical assays are more reproducible
        tech_reps = 2
        bio_reps = 2

    # Adjust for variant count
    if variant_count > 20:
        # Reduce replicates for high-throughput
        tech_reps = max(2, tech_reps - 1)
        bio_reps = max(2, bio_reps - 1)
    elif variant_count < 5:
        # Increase replicates for small studies
        tech_reps = min(5, tech_reps + 1)
        bio_reps = min(4, bio_reps + 1)

    return tech_reps, bio_reps


def _get_controls(work_package: WorkPackage) -> List[str]:
    """Generate list of control conditions."""
    controls = ["WT"]  # Wild-type control always included

    assay_module = work_package.assay_module

    if assay_module == "FUNCTIONAL":
        controls.extend(["no_substrate", "inhibitor_control"])
    elif assay_module == "TRAFFICKING":
        controls.extend(["ER_retention_control", "membrane_marker"])
    elif assay_module == "SPLICE_MINIGENE":
        controls.extend(["canonical_splice_control"])
    elif assay_module == "RNA_SEQ_REPORTER":
        controls.extend(["housekeeping_gene_control"])
    elif assay_module == "TRANSCRIPTIONAL_REPORTER":
        controls.extend(["empty_vector_control"])

    return controls


def create_design_summaries(work_packages: WorkPackages) -> DesignSummaries:
    """Create experimental designs for all work packages.

    Args:
        work_packages: List of WorkPackage objects

    Returns:
        List of DesignSummary objects
    """
    LOGGER.info("Creating experimental designs for %d work packages", len(work_packages))

    design_summaries = []

    for wp in work_packages:
        factors = _get_design_factors(wp)
        design_type = _get_design_type(factors, len(wp.variant_ids))
        tech_reps, bio_reps = _get_replicates(wp)
        controls = _get_controls(wp)

        design_summary = DesignSummary(
            wp_id=wp.wp_id,
            factors=factors,
            design_type=design_type,
            tech_reps=tech_reps,
            bio_reps=bio_reps,
            controls=controls,
        )
        design_summaries.append(design_summary)

    LOGGER.info("Created designs for %d work packages", len(design_summaries))
    return design_summaries


def save_design_summaries(
    design_summaries: DesignSummaries,
    output_dir: Path = CRO_DATA_DIR
) -> Path:
    """Save DesignSummaries to JSON file."""
    output_dir.mkdir(parents=True, exist_ok=True)

    json_path = output_dir / "design_summaries.json"
    with open(json_path, "w") as f:
        json.dump([{
            "wp_id": ds.wp_id,
            "factors": [{"name": f.name, "levels": f.levels} for f in ds.factors],
            "design_type": ds.design_type,
            "tech_reps": ds.tech_reps,
            "bio_reps": ds.bio_reps,
            "controls": ds.controls,
        } for ds in design_summaries], f, indent=2)

    LOGGER.info("Saved %d design summaries to %s", len(design_summaries), json_path)
    return json_path


def save_design_tables(
    work_packages: WorkPackages,
    design_summaries: DesignSummaries,
    output_dir: Path = CRO_DATA_DIR
) -> Path:
    """Save detailed design tables as CSV files for each work package with full Cartesian product matrices."""
    designs_dir = output_dir / "designs"
    designs_dir.mkdir(parents=True, exist_ok=True)

    # Create lookup for design summaries
    design_lookup = {ds.wp_id: ds for ds in design_summaries}

    for wp in work_packages:
        ds = design_lookup.get(wp.wp_id)
        if not ds:
            continue

        # Create base conditions (variants + controls)
        base_conditions = []

        # Add control conditions
        for control in ds.controls:
            base_conditions.append({
                "variant": control,
                "control_type": "experimental_control",
            })

        # Add variant conditions
        for variant_id in wp.variant_ids:
            base_conditions.append({
                "variant": variant_id,
                "control_type": "variant",
            })

        # Create factor combinations based on design type
        if ds.design_type == "one_factor":
            # For one-factor designs, use the first factor only
            if ds.factors:
                primary_factor = ds.factors[0]
                factor_levels = primary_factor.levels
            else:
                factor_levels = [None]  # No factors, just base conditions

            # Create all combinations of base conditions × factor levels
            design_combinations = list(itertools.product(base_conditions, factor_levels))

        elif ds.design_type == "full_factorial":
            # Full factorial: all combinations of all factors
            if ds.factors:
                factor_combinations = list(itertools.product(*[f.levels for f in ds.factors]))
            else:
                factor_combinations = [()]  # No factors

            # Create all combinations of base conditions × factor combinations
            design_combinations = list(itertools.product(base_conditions, factor_combinations))

        elif ds.design_type == "fractional":
            # Fractional design: sample a subset of the full factorial
            # For simplicity, sample 50% of the full factorial space or minimum 10 conditions
            if ds.factors:
                full_factorial = list(itertools.product(*[f.levels for f in ds.factors]))
                sample_size = max(10, len(full_factorial) // 2)
                # Use a seeded random sample for reproducibility
                random.seed(42)  # Fixed seed for reproducibility
                factor_combinations = random.sample(full_factorial, min(sample_size, len(full_factorial)))
                LOGGER.info(f"WP {wp.wp_id}: Fractional design sampled {len(factor_combinations)}/{len(full_factorial)} factor combinations ({sample_size} requested, {len(full_factorial)} available)")
            else:
                factor_combinations = [()]

            design_combinations = list(itertools.product(base_conditions, factor_combinations))

        else:
            # Fallback to one-factor behavior
            design_combinations = list(itertools.product(base_conditions, [None]))

        # Build design matrix
        design_data = []
        condition_counter = 0

        for base_condition, factor_values in design_combinations:
            condition_counter += 1
            condition_id = f"{wp.wp_id}_C{condition_counter:03d}"

            # Create row data
            row_data = {
                "wp_id": wp.wp_id,
                "condition_id": condition_id,
                "variant": base_condition["variant"],
                "control_type": base_condition["control_type"],
                "tech_reps": ds.tech_reps,
                "bio_reps": ds.bio_reps,
            }

            # Add factor columns
            if isinstance(factor_values, tuple):
                # Multiple factors
                for i, factor in enumerate(ds.factors):
                    if i < len(factor_values):
                        row_data[f"factor_{factor.name}"] = factor_values[i]
                    else:
                        row_data[f"factor_{factor.name}"] = None
            elif factor_values is not None and ds.factors:
                # Single factor
                row_data[f"factor_{ds.factors[0].name}"] = factor_values

            design_data.append(row_data)

        df = pd.DataFrame(design_data)

        # Note: Each row represents one experimental condition. The tech_reps and bio_reps columns
        # indicate how many technical and biological replicates to run for each condition.
        # CROs should expand these to individual replicate-level rows if needed.

        # Sanity check row counts for non-fractional designs
        # Each row represents one experimental condition (not including replicates)
        expected_conditions = len(base_conditions)
        if ds.design_type == "one_factor" and ds.factors:
            expected_conditions *= len(ds.factors[0].levels)
        elif ds.design_type == "full_factorial" and ds.factors:
            expected_conditions *= int(pd.Series([len(f.levels) for f in ds.factors]).product())

        if ds.design_type != "fractional" and len(df) != expected_conditions:
            LOGGER.warning(f"WP {wp.wp_id}: Design table has {len(df)} conditions, expected {expected_conditions} for {ds.design_type} design")
        else:
            LOGGER.info(f"WP {wp.wp_id}: Created design with {len(df)} experimental conditions")

        # Save to CSV
        csv_path = designs_dir / f"{wp.wp_id}_design.csv"
        df.to_csv(csv_path, index=False)
        LOGGER.info("Saved design table to %s", csv_path)

    return designs_dir


def main() -> None:
    """Run experimental design generation."""
    work_packages = load_work_packages()
    design_summaries = create_design_summaries(work_packages)

    save_design_summaries(design_summaries)
    save_design_tables(work_packages, design_summaries)

    LOGGER.info("Successfully created experimental designs for %d work packages", len(work_packages))


if __name__ == "__main__":
    main()
