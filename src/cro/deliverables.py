"""Stage 6: Deliverables Specification

Defines metrics, QC expectations, and deliverable manifests for each work package.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import List

try:
    # When run as module
    from .cro_types import AssayModuleId, DeliverableSpec, DeliverableSpecs, WorkPackages
    from .workpackages import load_work_packages
except ImportError:
    # When run as standalone script
    import sys
    from pathlib import Path
    # Add the src directory to sys.path
    src_dir = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(src_dir))
    from cro.cro_types import AssayModuleId, DeliverableSpec, DeliverableSpecs, WorkPackages
    from cro.workpackages import load_work_packages

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(name)s | %(levelname)s | %(message)s",
)
LOGGER = logging.getLogger(__name__)

CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]
CRO_DATA_DIR = CAMPAIGN_ROOT / "data_processed" / "cro"


def _get_primary_metrics(assay_module: AssayModuleId) -> List[str]:
    """Get primary metrics for an assay module."""
    metrics_map = {
        "DSF_SEC": [
            "melting_temperature_tm_celsius",
            "oligomeric_state_distribution",
            "thermal_stability_auc",
        ],
        "FUNCTIONAL": [
            "enzyme_activity_umol_min_mg",
            "kinetic_parameter_km_um",
            "kinetic_parameter_vmax",
            "inhibition_ic50_um",
        ],
        "TRAFFICKING": [
            "membrane_localization_score",
            "endoplasmic_reticulum_retention_score",
            "trafficking_efficiency_percent",
        ],
        "SPLICE_MINIGENE": [
            "splicing_efficiency_percent",
            "exon_inclusion_ratio",
            "cryptic_splice_site_usage",
        ],
        "RNA_SEQ_REPORTER": [
            "differential_expression_fold_change",
            "alternative_splicing_psi_score",
            "transcript_abundance_tpm",
        ],
        "TRANSCRIPTIONAL_REPORTER": [
            "promoter_activity_fold_induction",
            "transcription_factor_binding_affinity",
            "reporter_gene_expression_luminescence",
        ],
    }
    return metrics_map.get(assay_module, ["assay_result_score"])


def _get_raw_returns(assay_module: AssayModuleId) -> List[str]:
    """Get raw data return formats for an assay module."""
    returns_map = {
        "DSF_SEC": [
            "dsf_melting_curves_csv",
            "sec_chromatograms_pdf",
            "oligomer_analysis_report_pdf",
        ],
        "FUNCTIONAL": [
            "enzyme_kinetics_raw_data_csv",
            "activity_assay_plates_xlsx",
            "standard_curves_pdf",
        ],
        "TRAFFICKING": [
            "confocal_microscopy_images_tiff",
            "localization_quantification_csv",
            "image_analysis_report_pdf",
        ],
        "SPLICE_MINIGENE": [
            "splice_products_gel_images_tiff",
            "sequencing_validation_fastq",
            "minigene_construct_maps_pdf",
        ],
        "RNA_SEQ_REPORTER": [
            "raw_reads_fastq_gz",
            "alignment_bam_files",
            "gene_expression_counts_tsv",
        ],
        "TRANSCRIPTIONAL_REPORTER": [
            "luminescence_readings_csv",
            "dose_response_curves_pdf",
            "reporter_assay_plates_xlsx",
        ],
    }
    return returns_map.get(assay_module, ["raw_assay_data_csv"])


def _get_summary_columns(assay_module: AssayModuleId) -> List[str]:
    """Get summary data column specifications."""
    columns_map = {
        "DSF_SEC": [
            "variant_id",
            "condition",
            "replicate_id",
            "melting_temperature_celsius",
            "thermal_stability_score",
            "oligomeric_state",
            "quality_score",
        ],
        "FUNCTIONAL": [
            "variant_id",
            "condition",
            "replicate_id",
            "enzyme_activity",
            "substrate_concentration_um",
            "kinetic_parameter",
            "assay_quality_score",
        ],
        "TRAFFICKING": [
            "variant_id",
            "condition",
            "replicate_id",
            "compartment_localization",
            "trafficking_score",
            "membrane_insertion_efficiency",
            "image_quality_score",
        ],
        "SPLICE_MINIGENE": [
            "variant_id",
            "condition",
            "replicate_id",
            "splicing_efficiency",
            "exon_skipping_ratio",
            "cryptic_splicing_detected",
            "validation_confirmed",
        ],
        "RNA_SEQ_REPORTER": [
            "variant_id",
            "condition",
            "replicate_id",
            "gene_expression_fold_change",
            "splicing_event_psi",
            "transcript_abundance",
            "sequencing_depth",
        ],
        "TRANSCRIPTIONAL_REPORTER": [
            "variant_id",
            "condition",
            "replicate_id",
            "promoter_activity",
            "transcription_factor_binding",
            "reporter_expression_level",
            "assay_quality_score",
        ],
    }
    return columns_map.get(assay_module, ["variant_id", "condition", "replicate_id", "assay_result"])


def _get_qc_expectations(assay_module: AssayModuleId) -> List[str]:
    """Get QC expectations and criteria for an assay module."""
    qc_map = {
        "DSF_SEC": [
            "melting_temperature_std_dev_<_2c",
            "protein_concentration_0_5_2_mg_ml",
            "buffer_conditions_consistent",
            "thermal_ramp_rate_1c_min",
            "positive_controls_within_expected_range",
        ],
        "FUNCTIONAL": [
            "enzyme_activity_cv_<15_percent",
            "substrate_depletion_<20_percent",
            "standard_curve_r_squared_>_0_95",
            "blank_subtraction_performed",
            "positive_negative_controls_validated",
        ],
        "TRAFFICKING": [
            "cell_viability_>_80_percent",
            "transfection_efficiency_>_70_percent",
            "marker_colocalization_quantified",
            "image_resolution_adequate",
            "background_subtraction_applied",
        ],
        "SPLICE_MINIGENE": [
            "minigene_construct_integrity_confirmed",
            "transfection_efficiency_>_60_percent",
            "splice_products_resolved_by_gel",
            "sequencing_validation_performed",
            "canonical_splicing_control_valid",
        ],
        "RNA_SEQ_REPORTER": [
            "rna_integrity_rin_>_8",
            "sequencing_depth_>_10_million_reads",
            "mapping_rate_>_80_percent",
            "housekeeping_gene_expression_stable",
            "batch_effects_corrected",
        ],
        "TRANSCRIPTIONAL_REPORTER": [
            "cell_viability_>_85_percent",
            "transfection_efficiency_>_70_percent",
            "luminescence_signal_to_noise_>_10",
            "dose_response_curve_fitted",
            "background_subtraction_performed",
        ],
    }
    return qc_map.get(assay_module, ["assay_quality_controls_passed"])


def create_deliverable_specs(work_packages: WorkPackages) -> DeliverableSpecs:
    """Create deliverable specifications for all work packages.

    Args:
        work_packages: List of WorkPackage objects

    Returns:
        List of DeliverableSpec objects
    """
    LOGGER.info("Creating deliverable specifications for %d work packages", len(work_packages))

    deliverable_specs = []

    for wp in work_packages:
        primary_metrics = _get_primary_metrics(wp.assay_module)
        raw_returns = _get_raw_returns(wp.assay_module)
        summary_columns = _get_summary_columns(wp.assay_module)
        qc_expectations = _get_qc_expectations(wp.assay_module)

        spec = DeliverableSpec(
            wp_id=wp.wp_id,
            primary_metrics=primary_metrics,
            raw_returns=raw_returns,
            summary_columns=summary_columns,
            qc_expectations=qc_expectations,
        )
        deliverable_specs.append(spec)

    LOGGER.info("Created deliverable specs for %d work packages", len(deliverable_specs))
    return deliverable_specs


def save_deliverable_specs(
    deliverable_specs: DeliverableSpecs,
    output_dir: Path = CRO_DATA_DIR
) -> Path:
    """Save DeliverableSpecs to JSON file."""
    output_dir.mkdir(parents=True, exist_ok=True)

    json_path = output_dir / "deliverable_specs.json"
    with open(json_path, "w") as f:
        json.dump([{
            "wp_id": ds.wp_id,
            "primary_metrics": ds.primary_metrics,
            "raw_returns": ds.raw_returns,
            "summary_columns": ds.summary_columns,
            "qc_expectations": ds.qc_expectations,
        } for ds in deliverable_specs], f, indent=2)

    LOGGER.info("Saved %d deliverable specs to %s", len(deliverable_specs), json_path)
    return json_path


def main() -> None:
    """Run deliverable specification generation."""
    work_packages = load_work_packages()
    deliverable_specs = create_deliverable_specs(work_packages)
    save_deliverable_specs(deliverable_specs)

    LOGGER.info("Successfully created deliverable specifications for %d work packages", len(work_packages))


if __name__ == "__main__":
    main()
