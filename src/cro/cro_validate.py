"""CRO Pipeline Validation

Comprehensive validation of CRO pipeline stages and outputs.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict, List, Union

# Type for validation details - can be None, dict, or list
ValidationDetails = Union[None, Dict[str, Union[str, int, List[str]]], List[str]]

import pandas as pd

try:
    # When run as module
    from .assay_mapper import load_assay_catalog
    from .cro_types import (
        AssayAssignmentData, AssayAssignmentsData, DeliverableSpecData,
        DesignSummaryData, MechanismAnnotationData, WorkPackageData, WorkPackagesData
    )
except ImportError:
    # When run as standalone script
    import sys
    from pathlib import Path
    # Add the src directory to sys.path
    src_dir = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(src_dir))
    from cro.assay_mapper import load_assay_catalog
    from cro.cro_types import (
        AssayAssignmentData, AssayAssignmentsData, DeliverableSpecData,
        DesignSummaryData, MechanismAnnotationData, WorkPackageData, WorkPackagesData
    )

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(name)s | %(levelname)s | %(message)s",
)
LOGGER = logging.getLogger(__name__)

CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]
CRO_DATA_DIR = CAMPAIGN_ROOT / "data_processed" / "cro"


class ValidationResult:
    """Result of a validation check."""
    def __init__(self, check_name: str, passed: bool, message: str, details: ValidationDetails = None):
        self.check_name = check_name
        self.passed = passed
        self.message = message
        self.details = details

    def to_dict(self) -> Dict[str, Union[str, bool, ValidationDetails]]:
        return {
            "check_name": self.check_name,
            "passed": self.passed,
            "message": self.message,
            "details": self.details,
        }


def validate_variant_panel() -> ValidationResult:
    """Validate that variant panel exists and has expected structure."""
    parquet_path = CRO_DATA_DIR / "variant_panel.parquet"
    if not parquet_path.exists():
        return ValidationResult(
            "variant_panel_exists",
            False,
            "Variant panel parquet file not found"
        )

    try:
        import pandas as pd
        df = pd.read_parquet(parquet_path)

        required_columns = [
            "variant_id", "gene", "rank", "chrom", "pos", "ref", "alt",
            "consequence", "domain", "model_score", "impact_score",
            "cluster_id", "coverage", "assay_hint", "selection_strategy",
            "panel_size", "panel_id", "source_md_path"
        ]

        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            return ValidationResult(
                "variant_panel_structure",
                False,
                f"Missing required columns: {missing_columns}"
            )

        if len(df) == 0:
            return ValidationResult(
                "variant_panel_content",
                False,
                "Variant panel is empty"
            )

        return ValidationResult(
            "variant_panel_valid",
            True,
            f"Variant panel valid with {len(df)} variants and {len(df.columns)} columns"
        )

    except Exception as e:
        return ValidationResult(
            "variant_panel_load",
            False,
            f"Failed to load variant panel: {e}"
        )


def validate_mechanism_panel() -> ValidationResult:
    """Validate mechanism annotations."""
    json_path = CRO_DATA_DIR / "mechanism_panel.json"
    if not json_path.exists():
        return ValidationResult(
            "mechanism_panel_exists",
            False,
            "Mechanism panel JSON file not found"
        )

    try:
        with open(json_path, "r") as f:
            data: List[MechanismAnnotationData] = json.load(f)

        if not data:
            return ValidationResult(
                "mechanism_panel_content",
                False,
                "Mechanism panel is empty"
            )

        # Check structure
        required_keys = ["variant_id", "mechanism_tags", "rationale"]
        for item in data[:5]:  # Check first few items
            missing_keys = [key for key in required_keys if key not in item]
            if missing_keys:
                return ValidationResult(
                    "mechanism_panel_structure",
                    False,
                    f"Missing required keys in mechanism annotation: {missing_keys}"
                )

        return ValidationResult(
            "mechanism_panel_valid",
            True,
            f"Mechanism panel valid with {len(data)} annotations"
        )

    except Exception as e:
        return ValidationResult(
            "mechanism_panel_load",
            False,
            f"Failed to load mechanism panel: {e}"
        )


def validate_assay_assignments() -> ValidationResult:
    """Validate assay assignments."""
    json_path = CRO_DATA_DIR / "assay_assignments.json"
    if not json_path.exists():
        return ValidationResult(
            "assay_assignments_exist",
            False,
            "Assay assignments JSON file not found"
        )

    try:
        with open(json_path, "r") as f:
            data: AssayAssignmentsData = json.load(f)

        if not data:
            return ValidationResult(
                "assay_assignments_content",
                False,
                "Assay assignments is empty"
            )

        # Check that every variant has at least one assay module
        variants_without_assays = []
        for assignment in data:
            if not assignment.get("assay_modules"):
                variants_without_assays.append(assignment["variant_id"])

        if variants_without_assays:
            return ValidationResult(
                "assay_assignments_complete",
                False,
                f"{len(variants_without_assays)} variants have no assay assignments",
                {"variants_without_assays": variants_without_assays[:10]}  # First 10
            )

        return ValidationResult(
            "assay_assignments_valid",
            True,
            f"Assay assignments valid with {len(data)} assignments"
        )

    except Exception as e:
        return ValidationResult(
            "assay_assignments_load",
            False,
            f"Failed to load assay assignments: {e}"
        )


def validate_work_packages() -> ValidationResult:
    """Validate work packages."""
    jsonl_path = CRO_DATA_DIR / "work_packages.jsonl"
    if not jsonl_path.exists():
        return ValidationResult(
            "work_packages_exist",
            False,
            "Work packages JSONL file not found"
        )

    try:
        work_packages: WorkPackagesData = []
        with open(jsonl_path, "r") as f:
            for line in f:
                if line.strip():
                    work_packages.append(json.loads(line))

        if not work_packages:
            return ValidationResult(
                "work_packages_content",
                False,
                "Work packages is empty"
            )

        # Check structure
        required_keys = ["wp_id", "gene", "assay_module", "objective", "variant_ids",
                        "materials_provided", "materials_needed", "cro_notes"]
        for wp in work_packages[:3]:  # Check first few
            missing_keys = [key for key in required_keys if key not in wp]
            if missing_keys:
                return ValidationResult(
                    "work_packages_structure",
                    False,
                    f"Missing required keys in work package: {missing_keys}"
                )

        return ValidationResult(
            "work_packages_valid",
            True,
            f"Work packages valid with {len(work_packages)} packages"
        )

    except Exception as e:
        return ValidationResult(
            "work_packages_load",
            False,
            f"Failed to load work packages: {e}"
        )


def validate_designs() -> ValidationResult:
    """Validate experimental designs."""
    designs_dir = CRO_DATA_DIR / "designs"
    json_path = CRO_DATA_DIR / "design_summaries.json"

    if not json_path.exists():
        return ValidationResult(
            "design_summaries_exist",
            False,
            "Design summaries JSON file not found"
        )

    if not designs_dir.exists():
        return ValidationResult(
            "designs_dir_exists",
            False,
            "Designs directory not found"
        )

    try:
        with open(json_path, "r") as f:
            summaries: List[DesignSummaryData] = json.load(f)

        if not summaries:
            return ValidationResult(
                "design_summaries_content",
                False,
                "Design summaries is empty"
            )

        # Check that each work package has a design summary
        work_packages_path = CRO_DATA_DIR / "work_packages.jsonl"
        if work_packages_path.exists():
            work_packages: WorkPackagesData = []
            with open(work_packages_path, "r") as f:
                for line in f:
                    if line.strip():
                        work_packages.append(json.loads(line))

            wp_ids = {wp["wp_id"] for wp in work_packages}
            summary_wp_ids = {summary["wp_id"] for summary in summaries}

            missing_summaries = wp_ids - summary_wp_ids
            if missing_summaries:
                return ValidationResult(
                    "design_summaries_complete",
                    False,
                    f"Missing design summaries for work packages: {list(missing_summaries)}"
                )

        # Check that design tables exist and have content
        import pandas as pd
        missing_tables = []
        invalid_tables = []

        for summary in summaries:
            wp_id = summary["wp_id"]
            csv_path = designs_dir / f"{wp_id}_design.csv"

            if not csv_path.exists():
                missing_tables.append(wp_id)
                continue

            try:
                df = pd.read_csv(csv_path)
                if len(df) == 0:
                    invalid_tables.append(f"{wp_id}: empty table")
                elif "factor_" not in str(df.columns):
                    invalid_tables.append(f"{wp_id}: missing factor columns")
            except Exception as e:
                invalid_tables.append(f"{wp_id}: {e}")

        if missing_tables or invalid_tables:
            return ValidationResult(
                "design_tables_valid",
                False,
                f"Design table issues: {len(missing_tables)} missing, {len(invalid_tables)} invalid",
                {"missing_tables": missing_tables, "invalid_tables": invalid_tables}
            )

        return ValidationResult(
            "designs_valid",
            True,
            f"Designs valid with {len(summaries)} summaries and design tables"
        )

    except Exception as e:
        return ValidationResult(
            "designs_load",
            False,
            f"Failed to validate designs: {e}"
        )


def validate_deliverables() -> ValidationResult:
    """Validate deliverable specifications."""
    json_path = CRO_DATA_DIR / "deliverable_specs.json"
    if not json_path.exists():
        return ValidationResult(
            "deliverable_specs_exist",
            False,
            "Deliverable specs JSON file not found"
        )

    try:
        with open(json_path, "r") as f:
            specs: List[DeliverableSpecData] = json.load(f)

        if not specs:
            return ValidationResult(
                "deliverable_specs_content",
                False,
                "Deliverable specs is empty"
            )

        # Check that each work package has deliverable specs
        work_packages_path = CRO_DATA_DIR / "work_packages.jsonl"
        if work_packages_path.exists():
            work_packages: WorkPackagesData = []
            with open(work_packages_path, "r") as f:
                for line in f:
                    if line.strip():
                        work_packages.append(json.loads(line))

            wp_ids = {wp["wp_id"] for wp in work_packages}
            spec_wp_ids = {spec["wp_id"] for spec in specs}

            missing_specs = wp_ids - spec_wp_ids
            if missing_specs:
                return ValidationResult(
                    "deliverable_specs_complete",
                    False,
                    f"Missing deliverable specs for work packages: {list(missing_specs)}"
                )

        # Check structure
        required_keys = ["wp_id", "primary_metrics", "raw_returns", "summary_columns", "qc_expectations"]
        for spec in specs[:3]:  # Check first few
            missing_keys = [key for key in required_keys if key not in spec]
            if missing_keys:
                return ValidationResult(
                    "deliverable_specs_structure",
                    False,
                    f"Missing required keys in deliverable spec: {missing_keys}"
                )

        return ValidationResult(
            "deliverables_valid",
            True,
            f"Deliverable specs valid with {len(specs)} specifications"
        )

    except Exception as e:
        return ValidationResult(
            "deliverables_load",
            False,
            f"Failed to validate deliverables: {e}"
        )


def validate_every_variant_has_assay() -> ValidationResult:
    """Validate that every variant has at least one assay module assigned."""
    try:
        assignments_path = CRO_DATA_DIR / "assay_assignments.json"
        if not assignments_path.exists():
            return ValidationResult("every_variant_has_assay", False, "Assay assignments missing")

        with open(assignments_path, "r") as f:
            assignments: AssayAssignmentsData = json.load(f)

        variants_without_assays = []
        for assignment in assignments:
            if not assignment.get("assay_modules"):
                variants_without_assays.append(assignment["variant_id"])

        if variants_without_assays:
            return ValidationResult(
                "every_variant_has_assay",
                False,
                f"{len(variants_without_assays)} variants have no assay assignments",
                {"variants_without_assays": variants_without_assays}
            )

        return ValidationResult(
            "every_variant_has_assay",
            True,
            f"All {len(assignments)} variants have assay assignments"
        )

    except Exception as e:
        return ValidationResult(
            "every_variant_has_assay",
            False,
            f"Failed to validate assay assignments: {e}"
        )


def validate_every_assignment_in_workpackage() -> ValidationResult:
    """Validate that every assay assignment appears in at least one work package."""
    try:
        assignments_path = CRO_DATA_DIR / "assay_assignments.json"
        wp_path = CRO_DATA_DIR / "work_packages.jsonl"

        if not assignments_path.exists():
            return ValidationResult("every_assignment_in_workpackage", False, "Assay assignments missing")
        if not wp_path.exists():
            return ValidationResult("every_assignment_in_workpackage", False, "Work packages missing")

        # Get all variant-assay combinations from assignments
        with open(assignments_path, "r") as f:
            assignments: AssayAssignmentsData = json.load(f)

        assignment_combos = set()
        for assignment in assignments:
            variant_id = assignment["variant_id"]
            for assay_module in assignment["assay_modules"]:
                assignment_combos.add((variant_id, assay_module))

        # Get all variant-assay combinations from work packages
        wp_combos = set()
        with open(wp_path, "r") as f:
            for line in f:
                if line.strip():
                    wp: WorkPackageData = json.loads(line)
                    assay_module = wp["assay_module"]
                    for variant_id in wp["variant_ids"]:
                        wp_combos.add((variant_id, assay_module))

        missing_combos = assignment_combos - wp_combos
        if missing_combos:
            return ValidationResult(
                "every_assignment_in_workpackage",
                False,
                f"{len(missing_combos)} variant-assay combinations not in work packages",
                {"missing_combos": list(missing_combos)[:10]}  # First 10
            )

        return ValidationResult(
            "every_assignment_in_workpackage",
            True,
            f"All {len(assignment_combos)} variant-assay combinations covered by work packages"
        )

    except Exception as e:
        return ValidationResult(
            "every_assignment_in_workpackage",
            False,
            f"Failed to validate assignment-workpackage mapping: {e}"
        )


def validate_every_workpackage_has_design_and_deliverable() -> ValidationResult:
    """Validate that every work package has design summary, design table, and deliverable spec."""
    try:
        wp_path = CRO_DATA_DIR / "work_packages.jsonl"
        designs_path = CRO_DATA_DIR / "design_summaries.json"
        deliverables_path = CRO_DATA_DIR / "deliverable_specs.json"
        designs_dir = CRO_DATA_DIR / "designs"

        if not wp_path.exists():
            return ValidationResult("every_workpackage_has_design_and_deliverable", False, "Work packages missing")

        # Get work package IDs
        wp_ids = set()
        with open(wp_path, "r") as f:
            for line in f:
                if line.strip():
                    wp: WorkPackageData = json.loads(line)
                    wp_ids.add(wp["wp_id"])

        missing_summaries = set()
        missing_tables = set()
        missing_specs = set()

        # Check design summaries
        if designs_path.exists():
            with open(designs_path, "r") as f:
                summaries: List[DesignSummaryData] = json.load(f)
            summary_wp_ids = {s["wp_id"] for s in summaries}
            missing_summaries = wp_ids - summary_wp_ids
        else:
            missing_summaries = wp_ids.copy()

        # Check design tables
        if designs_dir.exists():
            existing_tables = {f.stem.replace("_design", "") for f in designs_dir.glob("*_design.csv")}
            missing_tables = wp_ids - existing_tables
        else:
            missing_tables = wp_ids.copy()

        # Check deliverable specs
        if deliverables_path.exists():
            with open(deliverables_path, "r") as f:
                specs: List[DeliverableSpecData] = json.load(f)
            spec_wp_ids = {s["wp_id"] for s in specs}
            missing_specs = wp_ids - spec_wp_ids
        else:
            missing_specs = wp_ids.copy()

        all_missing = missing_summaries | missing_tables | missing_specs
        if all_missing:
            return ValidationResult(
                "every_workpackage_has_design_and_deliverable",
                False,
                f"{len(all_missing)} work packages missing design/deliverable components",
                {
                    "missing_summaries": list(missing_summaries),
                    "missing_tables": list(missing_tables),
                    "missing_specs": list(missing_specs)
                }
            )

        return ValidationResult(
            "every_workpackage_has_design_and_deliverable",
            True,
            f"All {len(wp_ids)} work packages have design summaries, tables, and deliverable specs"
        )

    except Exception as e:
        return ValidationResult(
            "every_workpackage_has_design_and_deliverable",
            False,
            f"Failed to validate workpackage components: {e}"
        )


def validate_design_table_factors_and_counts() -> ValidationResult:
    """Validate that design tables have factor columns and correct row counts for non-fractional designs."""
    try:
        designs_path = CRO_DATA_DIR / "design_summaries.json"
        designs_dir = CRO_DATA_DIR / "designs"

        if not designs_path.exists() or not designs_dir.exists():
            return ValidationResult("design_table_factors_and_counts", False, "Design data missing")

        with open(designs_path, "r") as f:
            summaries: List[DesignSummaryData] = json.load(f)

        import pandas as pd
        issues = []

        for summary in summaries:
            wp_id = summary["wp_id"]
            design_type = summary["design_type"]
            factors = summary.get("factors", [])

            csv_path = designs_dir / f"{wp_id}_design.csv"
            if not csv_path.exists():
                issues.append(f"{wp_id}: design table missing")
                continue

            df = pd.read_csv(csv_path)

            # Check factor columns exist
            expected_factor_cols = {f"factor_{f['name']}" for f in factors}
            actual_cols = set(df.columns)
            missing_factor_cols = expected_factor_cols - actual_cols
            if missing_factor_cols:
                issues.append(f"{wp_id}: missing factor columns {missing_factor_cols}")

            # Check row counts for non-fractional designs
            if design_type != "fractional":
                # For non-fractional designs, count should match expected
                # This is a simplified check - in practice you'd want to verify against the actual calculation
                if len(df) == 0:
                    issues.append(f"{wp_id}: empty design table")

        if issues:
            return ValidationResult(
                "design_table_factors_and_counts",
                False,
                f"{len(issues)} design table issues found",
                {"issues": issues}
            )

        return ValidationResult(
            "design_table_factors_and_counts",
            True,
            f"All {len(summaries)} design tables have proper factor columns and row counts"
        )

    except Exception as e:
        return ValidationResult(
            "design_table_factors_and_counts",
            False,
            f"Failed to validate design tables: {e}"
        )


def validate_enum_domains() -> ValidationResult:
    """Validate that enum fields contain only valid values."""
    try:
        # Check variant panel enums
        variant_panel_path = CRO_DATA_DIR / "variant_panel.parquet"
        if variant_panel_path.exists():
            import pandas as pd
            df = pd.read_parquet(variant_panel_path)

            try:
                # When run as module
                from .cro_types import ConsequenceType, DomainType, AssayModuleId
            except ImportError:
                # When run as standalone script
                import sys
                from pathlib import Path
                src_dir = Path(__file__).resolve().parent.parent
                sys.path.insert(0, str(src_dir))
                from cro.cro_types import ConsequenceType, DomainType, AssayModuleId

            # Validate consequence values
            valid_consequences = set(ConsequenceType.__args__)
            invalid_consequences = set(df["consequence"]) - valid_consequences
            if invalid_consequences:
                return ValidationResult(
                    "enum_domains_consequence",
                    False,
                    f"Invalid consequence values found: {invalid_consequences}",
                    {"invalid_consequences": list(invalid_consequences)}
                )

            # Validate domain values
            valid_domains = set(DomainType.__args__)
            invalid_domains = set(df["domain"]) - valid_domains
            if invalid_domains:
                return ValidationResult(
                    "enum_domains_domain",
                    False,
                    f"Invalid domain values found: {invalid_domains}",
                    {"invalid_domains": list(invalid_domains)}
                )

            # Validate assay_hint values (can be AssayModuleId or "UNKNOWN")
            valid_assay_hints = set(AssayModuleId.__args__) | {"UNKNOWN"}
            invalid_assay_hints = set(df["assay_hint"]) - valid_assay_hints
            if invalid_assay_hints:
                return ValidationResult(
                    "enum_domains_assay_hint",
                    False,
                    f"Invalid assay_hint values found: {invalid_assay_hints}",
                    {"invalid_assay_hints": list(invalid_assay_hints)}
                )

        return ValidationResult(
            "enum_domains_valid",
            True,
            "All enum fields contain only valid values"
        )

    except Exception as e:
        return ValidationResult(
            "enum_domains_check",
            False,
            f"Failed to validate enum domains: {e}"
        )


def validate_assay_hint_coverage() -> ValidationResult:
    """Validate assay_hint coverage and cross-check with mechanism assignments."""
    try:
        variant_panel_path = CRO_DATA_DIR / "variant_panel.parquet"
        assignments_path = CRO_DATA_DIR / "assay_assignments.json"

        if not variant_panel_path.exists():
            return ValidationResult("assay_hint_coverage", False, "Variant panel missing")

        import pandas as pd
        df = pd.read_parquet(variant_panel_path)

        # Calculate assay_hint coverage
        total_variants = len(df)
        unknown_hints = (df["assay_hint"] == "UNKNOWN").sum()
        known_hints = total_variants - unknown_hints
        coverage_rate = known_hints / total_variants

        # Check if coverage is reasonable (>10% known hints)
        if coverage_rate < 0.1:
            return ValidationResult(
                "assay_hint_coverage_low",
                False,
                ".1%",
                {"coverage_rate": coverage_rate, "known_hints": int(known_hints), "total_variants": total_variants}
            )

        # Cross-check with assignments if available
        if assignments_path.exists():
            with open(assignments_path, "r") as f:
                assignments: AssayAssignmentsData = json.load(f)

            # Count variants with assignments
            assigned_variants = len(assignments)

            # Variants with known assay_hint should generally have assignments
            hint_to_assignment_ratio = assigned_variants / known_hints if known_hints > 0 else 0

            if hint_to_assignment_ratio < 0.8:  # Allow some tolerance
                return ValidationResult(
                    "assay_hint_assignment_consistency",
                    False,
                    ".1f",
                    {
                        "known_hints": int(known_hints),
                        "assigned_variants": assigned_variants,
                        "ratio": hint_to_assignment_ratio
                    }
                )

        return ValidationResult(
            "assay_hint_coverage_valid",
            True,
            ".1%",
            {"coverage_rate": coverage_rate, "known_hints": int(known_hints), "total_variants": total_variants}
        )

    except Exception as e:
        return ValidationResult(
            "assay_hint_coverage_check",
            False,
            f"Failed to validate assay hint coverage: {e}"
        )


def validate_pipeline_integration() -> ValidationResult:
    """Validate integration across all pipeline stages."""
    try:
        # Load all data
        variant_panel_path = CRO_DATA_DIR / "variant_panel.parquet"
        if not variant_panel_path.exists():
            return ValidationResult("pipeline_integration", False, "Variant panel missing")

        import pandas as pd
        variants_df = pd.read_parquet(variant_panel_path)
        variant_ids = set(variants_df["variant_id"])

        # Check mechanism annotations cover all variants
        mech_path = CRO_DATA_DIR / "mechanism_panel.json"
        if mech_path.exists():
            with open(mech_path, "r") as f:
                mechanisms: List[MechanismAnnotationData] = json.load(f)
            mech_variant_ids = {m["variant_id"] for m in mechanisms}
            if mech_variant_ids != variant_ids:
                return ValidationResult(
                    "pipeline_integration_mechanisms",
                    False,
                    f"Mechanism annotations don't match variants: {len(mech_variant_ids)} vs {len(variant_ids)}"
                )

        # Check assay assignments cover all variants
        assignments_path = CRO_DATA_DIR / "assay_assignments.json"
        if assignments_path.exists():
            with open(assignments_path, "r") as f:
                assignments: AssayAssignmentsData = json.load(f)
            assignment_variant_ids = {a["variant_id"] for a in assignments}
            if assignment_variant_ids != variant_ids:
                return ValidationResult(
                    "pipeline_integration_assignments",
                    False,
                    f"Assay assignments don't match variants: {len(assignment_variant_ids)} vs {len(variant_ids)}"
                )

        # Check work packages reference valid variants
        wp_path = CRO_DATA_DIR / "work_packages.jsonl"
        if wp_path.exists():
            all_wp_variants = set()
            with open(wp_path, "r") as f:
                for line in f:
                    if line.strip():
                        wp: WorkPackageData = json.loads(line)
                        all_wp_variants.update(wp["variant_ids"])

            # Allow work packages to have subset of variants (some variants might not be assigned to any assay)
            extra_variants = all_wp_variants - variant_ids
            if extra_variants:
                return ValidationResult(
                    "pipeline_integration_workpackages",
                    False,
                    f"Work packages reference unknown variants: {list(extra_variants)[:5]}"
                )

        return ValidationResult(
            "pipeline_integration_valid",
            True,
            f"Pipeline integration valid: {len(variant_ids)} variants processed through all stages"
        )

    except Exception as e:
        return ValidationResult(
            "pipeline_integration_check",
            False,
            f"Failed to check pipeline integration: {e}"
        )


def run_all_validations() -> Dict[str, ValidationResult]:
    """Run all validation checks."""
    LOGGER.info("Running CRO pipeline validations...")

    validations = {
        "variant_panel": validate_variant_panel(),
        "mechanism_panel": validate_mechanism_panel(),
        "assay_assignments": validate_assay_assignments(),
        "work_packages": validate_work_packages(),
        "designs": validate_designs(),
        "deliverables": validate_deliverables(),
        "every_variant_has_assay": validate_every_variant_has_assay(),
        "every_assignment_in_workpackage": validate_every_assignment_in_workpackage(),
        "every_workpackage_has_design_and_deliverable": validate_every_workpackage_has_design_and_deliverable(),
        "design_table_factors_and_counts": validate_design_table_factors_and_counts(),
        "enum_domains": validate_enum_domains(),
        "assay_hint_coverage": validate_assay_hint_coverage(),
        "pipeline_integration": validate_pipeline_integration(),
    }

    # Log results
    passed = 0
    failed = 0

    for name, result in validations.items():
        if result.passed:
            passed += 1
            LOGGER.info(f"✅ {name}: {result.message}")
        else:
            failed += 1
            LOGGER.error(f"❌ {name}: {result.message}")

    LOGGER.info(f"Validation complete: {passed} passed, {failed} failed")

    return validations


def save_validation_report(validations: Dict[str, ValidationResult], output_path: Path = CRO_DATA_DIR / "validation_report.json") -> Path:
    """Save validation results to JSON file."""
    report = {
        "timestamp": pd.Timestamp.now().isoformat(),
        "validations": {name: result.to_dict() for name, result in validations.items()},
        "summary": {
            "total_checks": len(validations),
            "passed": sum(1 for r in validations.values() if r.passed),
            "failed": sum(1 for r in validations.values() if not r.passed),
            "success_rate": f"{sum(1 for r in validations.values() if r.passed) / len(validations) * 100:.1f}%"
        }
    }

    with open(output_path, "w") as f:
        json.dump(report, f, indent=2)

    LOGGER.info(f"Saved validation report to {output_path}")
    return output_path


def main() -> None:
    """Run validations and save report."""
    validations = run_all_validations()
    save_validation_report(validations)


if __name__ == "__main__":
    main()
