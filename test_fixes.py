#!/usr/bin/env python3
"""
Test script to validate the implemented fixes according to the junior-friendly checklist.
"""

import pandas as pd
import numpy as np
import sys
from pathlib import Path

# Add the project root to Python path
sys.path.insert(0, str(Path(__file__).resolve().parents[0]))

def test_feature_engineering_fixes():
    """Test the feature engineering fixes according to the validation checklist."""
    print("🧪 Testing Feature Engineering Fixes")
    print("=" * 50)

    # Load the final features file
    features_path = Path("data_processed/features/variants_features_raw.parquet")
    if not features_path.exists():
        print("❌ variants_features_raw.parquet not found. Run feature engineering first.")
        return False

    df = pd.read_parquet(features_path)
    print(f"✅ Loaded {len(df)} variants with {len(df.columns)} columns")

    # Test 1: Check for minimal row loss (allow small loss due to deduplication)
    annotated_path = Path("data_processed/annotations/abca4_vus_annotated.parquet")
    if annotated_path.exists():
        df_annotated = pd.read_parquet(annotated_path)
        row_loss = len(df_annotated) - len(df)
        if row_loss <= 10:  # Allow small loss due to deduplication
            print(f"✅ (a) Acceptable row loss: {len(df_annotated)} → {len(df)} ({row_loss} rows lost, likely deduplication)")
        else:
            print(f"❌ (a) Excessive row loss: {len(df_annotated)} → {len(df)} ({row_loss} rows lost)")
            return False
    else:
        print("⚠️ Cannot check row loss - annotated file not found")

    # Test 2: Check key feature non-null rates (realistic thresholds)
    key_columns_checks = {
        'alphamissense_score': ('AlphaMissense', 40),  # Only works on missense variants
        'spliceai_max_score': ('SpliceAI', 95),        # Works on most variants
        'phylop_score': ('conservation', 95),          # Works on most variants
        'gnomad_max_af': ('AF (max)', 95),             # Population frequency data
        'domain': ('domain', 95),                      # Domain assignment
    }

    # Check for any available AF columns
    af_columns = [c for c in df.columns if 'af' in c.lower() and c != 'af_penalty']
    if af_columns:
        key_columns_checks[af_columns[0]] = ('AF (available)', 95)

    all_pass = True
    for col, (desc, min_pct) in key_columns_checks.items():
        if col in df.columns:
            non_null_pct = (df[col].notna().sum() / len(df)) * 100
            if non_null_pct >= min_pct:
                print(f"✅ {desc}: {non_null_pct:.1f}% non-null")
            else:
                print(f"❌ {desc}: {non_null_pct:.1f}% non-null (<{min_pct}% threshold)")
                all_pass = False
        else:
            print(f"❌ {desc}: Column '{col}' not found")
            all_pass = False

    # Test 3: Check AF penalty is non-zero for common variants
    if 'af_penalty' in df.columns and 'gnomad_max_af' in df.columns:
        common_variants = df[df['gnomad_max_af'] > 0.01]  # Common variants
        if not common_variants.empty:
            non_zero_penalty = (common_variants['af_penalty'] != 0).sum()
            if non_zero_penalty > 0:
                print(f"✅ (c) AF penalty non-zero for {non_zero_penalty}/{len(common_variants)} common variants")
            else:
                print("❌ (c) AF penalty is zero for all common variants")
                all_pass = False
        else:
            print("⚠️ No common variants found for AF penalty test")
    else:
        print("❌ Missing af_penalty or gnomad_max_af columns")
        all_pass = False

    # Test 4: Check domain_flag is properly set (should be 1 for all since all variants get domain assignment)
    if 'domain_flag' in df.columns:
        ones_count = (df['domain_flag'] == 1).sum()
        zeros_count = (df['domain_flag'] == 0).sum()
        if ones_count == len(df) and zeros_count == 0:
            print(f"✅ (d) domain_flag correctly set: {ones_count} ones (all variants have domain assignment)")
        else:
            print(f"❌ (d) domain_flag incorrect: {zeros_count} zeros, {ones_count} ones")
            all_pass = False
    else:
        print("❌ domain_flag column not found")
        all_pass = False

    return all_pass

def test_clustering_fixes():
    """Test the clustering fixes."""
    print("\n🧬 Testing Clustering Fixes")
    print("=" * 50)

    # Check domain computation
    domain_path = Path("data_processed/features/variants_with_domains.parquet")
    if not domain_path.exists():
        print("❌ variants_with_domains.parquet not found")
        return False

    df_domains = pd.read_parquet(domain_path)
    domain_counts = df_domains['domain'].value_counts()

    # Check TMD split
    tmd_domains = [d for d in domain_counts.index if 'TMD' in d]
    if 'TMD1' in tmd_domains and 'TMD2' in tmd_domains:
        print("✅ TMD split into TMD1/TMD2")
    else:
        print(f"❌ TMD not properly split: found {tmd_domains}")
        return False

    # Check τⱼ recalibration would work (can't fully test without running notebook)
    print("✅ Domain computation includes TMD1/TMD2 split")

    return True

def test_llm_assay_drafts():
    """Test LLM assay draft generation components."""
    print("\n🤖 Testing LLM Assay Drafts")
    print("=" * 50)

    import os
    import json
    import tempfile
    from unittest.mock import patch, MagicMock
    from pathlib import Path

    # Import the modules we need to test
    try:
        from src.config import REQUIRED_VARIANT_COLUMNS, validate_dataframe, load_parquet_safely
        from src.reporting.generate_assay_drafts import (
            validate_selected_variants_panel, prepare_variant_for_prompt
        )
        from src.reporting.llm_client import get_prompt_hash, validate_assay_response
    except ImportError as e:
        print(f"❌ Cannot import LLM modules: {e}")
        return False

    all_pass = True

    # Test 1: Data contract validation - required columns enforcement
    print("Testing data contract validation...")
    test_df = pd.DataFrame({
        'variant_id': ['test1', 'test2'],
        'gene': ['ABCA4', 'ABCA4'],
        'vep_consequence': ['missense_variant', 'missense_variant'],
        'cluster_id': ['TMD1_missense_variant', 'TMD2_missense_variant'],
        'impact_score': [0.8, 0.7],
        'gnomad_max_af': [0.001, 0.002]
    })

    try:
        validate_selected_variants_panel(test_df)
        print("✅ Required columns validation passed")
    except Exception as e:
        print(f"❌ Required columns validation failed: {e}")
        all_pass = False

    # Test missing required column
    bad_df = test_df.drop('vep_consequence', axis=1)
    try:
        validate_selected_variants_panel(bad_df)
        print("❌ Should have failed on missing required column")
        all_pass = False
    except ValueError:
        print("✅ Correctly rejected missing required column")

    # Test 2: Variant preparation for prompt
    print("Testing variant preparation...")
    row = test_df.iloc[0]
    variant_dict = prepare_variant_for_prompt(row)

    expected_keys = [
        'variant_id', 'gene', 'protein_change', 'consequence',
        'domain_cluster_id', 'mechanism', 'impact_score', 'gnomad_max_af'
    ]

    if all(key in variant_dict for key in expected_keys):
        print("✅ Variant preparation includes all required fields")
    else:
        print("❌ Variant preparation missing fields")
        all_pass = False

    # Test 3: Prompt hash generation (deterministic)
    print("Testing prompt hash generation...")
    test_prompt = "Test prompt for assay generation"
    hash1 = get_prompt_hash(test_prompt)
    hash2 = get_prompt_hash(test_prompt)

    if hash1 == hash2 and len(hash1) == 16:
        print("✅ Prompt hash generation is deterministic and correct length")
    else:
        print("❌ Prompt hash generation failed")
        all_pass = False

    # Test 4: Assay response validation
    print("Testing assay response validation...")
    valid_response = """
    Assay Type: minigene

    Cell Line: HEK293T

    Construct Design: Full length ABCA4 in pcDNA3.1

    Readout: RT-PCR quantification

    Controls: WT plasmid, empty vector

    Expected WT vs Mutant: Reduced splicing for mutant

    Effort: M

    Budget Note: Standard molecular biology costs

    Protocol: 1. Transfect cells with plasmids at 70% confluence. 2. Incubate for 48 hours. 3. Extract total RNA using Trizol reagent. 4. Perform RT-PCR with primers targeting the exon junction. 5. Quantify band intensities using gel densitometry. This is a complete protocol that should work for testing splicing defects in ABCA4 variants.
    """

    if validate_assay_response(valid_response):
        print("✅ Valid assay response correctly validated")
    else:
        print("❌ Valid assay response rejected")
        all_pass = False

    invalid_response = "This is not a proper assay response."
    if not validate_assay_response(invalid_response):
        print("✅ Invalid assay response correctly rejected")
    else:
        print("❌ Invalid assay response incorrectly accepted")
        all_pass = False

    # Test 5: Static prompt template rendering
    print("Testing prompt template rendering...")
    try:
        from jinja2 import Environment, FileSystemLoader
        from src.config import CAMPAIGN_ROOT

        template_dir = CAMPAIGN_ROOT / "src" / "reporting" / "templates"
        env = Environment(loader=FileSystemLoader(template_dir))
        template = env.get_template("assay_prompt.md.jinja")

        rendered = template.render(**variant_dict)
        word_count = len(rendered.split())

        if word_count < 300:  # Should be well under limit
            print(f"✅ Template renders correctly ({word_count} words)")
        else:
            print(f"❌ Template too long ({word_count} words)")
            all_pass = False

        # Check that variables were substituted
        if variant_dict['variant_id'] in rendered:
            print("✅ Template variables correctly substituted")
        else:
            print("❌ Template variables not substituted")
            all_pass = False

    except Exception as e:
        print(f"❌ Template rendering failed: {e}")
        all_pass = False

    return all_pass

def main():
    """Run all tests."""
    print("🔬 ABCA4 Pipeline Fix Validation")
    print("=" * 60)

    success = True

    # Test feature engineering fixes
    if not test_feature_engineering_fixes():
        success = False

    # Test clustering fixes
    if not test_clustering_fixes():
        success = False

    # Test LLM assay drafts
    if not test_llm_assay_drafts():
        success = False

    print("\n" + "=" * 60)
    if success:
        print("🎉 All validation checks passed!")
        print("\n📋 Next steps:")
        print("1. Run notebooks/02_feature_engineering.py fully")
        print("2. Inspect QC plots for AUROC ≥0.8 and P/LP right-shift")
        print("3. Run notebooks/03_optimization_dashboard.py with K=30, λ=0.6")
        print("4. Verify comparison panel shows similar total impact (<5% drop)")
        print("5. Check cluster coverage improvements")
    else:
        print("❌ Some validation checks failed. Review and fix issues above.")
        sys.exit(1)

if __name__ == "__main__":
    main()
