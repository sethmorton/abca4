#!/usr/bin/env python3
"""
Master script to run the full MAVE ingestion, normalization, and evaluation pipeline.

Usage:
    python src/mave/run_mave_pipeline.py --phase all
    python src/mave/run_mave_pipeline.py --phase ingest
    python src/mave/run_mave_pipeline.py --check
"""

import logging
import sys
import argparse
from pathlib import Path
from typing import Optional, Dict

import pandas as pd

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Add project root to path
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.mave.pipeline import ingest, normalize, features, strategies
from src.mave.evaluation import eval, sanity


def run_phase_ingest(
    datasets_config: Path = Path("config/mave_datasets.yaml"),
    output_base: Path = Path("data_processed/mave"),
) -> dict:
    """Run Phase B: Ingestion."""
    logger.info("\n" + "=" * 70)
    logger.info("PHASE B: INGEST")
    logger.info("=" * 70)
    
    results = ingest.ingest_all_datasets(datasets_config, output_base)
    
    logger.info(f"Ingestion complete: {len(results)} datasets ingested")
    return results


def run_phase_normalize(
    datasets_config: Path = Path("config/mave_datasets.yaml"),
    base_dir: Path = Path("data_processed/mave"),
) -> dict:
    """Run Phase C: Normalization and hit definition."""
    logger.info("\n" + "=" * 70)
    logger.info("PHASE C: NORMALIZE")
    logger.info("=" * 70)
    
    results = normalize.normalize_all_datasets(datasets_config, base_dir)
    
    logger.info(f"Normalization complete: {len(results)} datasets normalized")
    return results


def run_phase_features(
    base_dir: Path = Path("data_processed/mave"),
    datasets_config: Path = Path("config/mave_datasets.yaml"),
    gene_config: Optional[Dict] = None,
) -> dict:
    """Run Phase D: Feature pipeline integration."""
    logger.info("\n" + "=" * 70)
    logger.info("PHASE D: FEATURES")
    logger.info("=" * 70)
    
    import yaml
    from src.mave.pipeline.features import create_dms_full
    
    if not datasets_config.exists():
        logger.warning(f"Config not found: {datasets_config}")
        return {}
    
    try:
        with open(datasets_config) as f:
            config = yaml.safe_load(f)
    except Exception as e:
        logger.error(f"Failed to load config: {e}")
        return {}
    
    datasets = config.get("datasets", [])
    if not datasets:
        logger.warning("No datasets in config")
        return {}
    
    results = {}
    for dataset_cfg in datasets:
        dataset_id = dataset_cfg.get("dataset_id")
        if not dataset_id:
            continue
        
        dataset_dir = base_dir / dataset_id
        gene = dataset_cfg.get("gene")
        
        logger.info(f"\nProcessing {dataset_id}...")
        df = create_dms_full(dataset_id, dataset_dir, gene, gene_config)
        
        if df is not None:
            results[dataset_id] = df
    
    return results


def run_phase_evaluation(
    base_dir: Path = Path("data_processed/mave"),
    k_values: list = None,
    results_dir: Path = Path("results/mave"),
) -> dict:
    """Run Phase E: Selection strategies and evaluation."""
    if k_values is None:
        k_values = [10]
    
    logger.info("\n" + "=" * 70)
    logger.info("PHASE E: EVALUATION")
    logger.info("=" * 70)
    
    from src.mave.pipeline.strategies import run_all_strategies
    from src.mave.evaluation.eval import evaluate_dataset, save_metrics, sanity_check_metrics
    
    results_dir.mkdir(parents=True, exist_ok=True)
    all_results = {}
    
    # Find all dataset directories
    dataset_dirs = [d for d in base_dir.iterdir() if d.is_dir()]
    
    for dataset_dir in sorted(dataset_dirs):
        dataset_id = dataset_dir.name
        full_parquet = dataset_dir / f"{dataset_id}_full.parquet"
        
        if not full_parquet.exists():
            logger.warning(f"Skipping {dataset_id}: no full parquet")
            continue
        
        try:
            df_full = pd.read_parquet(full_parquet)
            logger.info(f"\nEvaluating {dataset_id}: {len(df_full)} variants, {df_full['is_hit'].sum()} hits")
        except Exception as e:
            logger.error(f"Failed to load {full_parquet}: {e}")
            continue
        
        for k in k_values:
            if k > len(df_full):
                logger.warning(f"k={k} > n={len(df_full)}, skipping")
                continue
            
            logger.info(f"  Running strategies for k={k}...")
            
            # Import strategies and Strand optimizer
            from src.mave.pipeline import strategies
            from src.reward.optimization import VariantOptimizer
            
            # Run all strategies including Strand
            strategy_results = {}
            
            # Baselines
            strategy_results["random"] = strategies.select_random(df_full, k)
            strategy_results["top_model_score"] = strategies.select_top_model_score(df_full, k, score_col="impact_score")
            strategy_results["top_conservation"] = strategies.select_top_conservation(df_full, k, conservation_col="conservation")
            strategy_results["oracle_functional"] = strategies.select_top_functional_score(df_full, k, score_col="functional_score")
            
            # Your Strand selection algorithm (REAL, not synthetic)
            try:
                strategy_results["strand"] = strategies.select_strand(
                    df_full,
                    k,
                    config={},
                    select_greedy_func=VariantOptimizer.select_greedy,
                    lambda_penalty=0.6
                )
                logger.info(f"  strand: OK")
            except Exception as e:
                logger.error(f"  strand: FAILED - {str(e)[:60]}")
                raise  # Fail fast - no fallback
            
            # Evaluate
            metrics_df = evaluate_dataset(dataset_id, df_full, strategy_results, [k])
            
            # Sanity checks
            issues = sanity_check_metrics(metrics_df, dataset_id)
            if issues:
                logger.warning(f"  Sanity check issues for {dataset_id}:")
                for issue in issues:
                    logger.warning(f"    - {issue}")
            
            # Save
            save_metrics(metrics_df, f"{dataset_id}_k{k}", results_dir)
            all_results[f"{dataset_id}_k{k}"] = metrics_df
    
    return all_results


def run_sanity_checks(
    base_dir: Path = Path("data_processed/mave"),
    datasets_config: Path = Path("config/mave_datasets.yaml"),
) -> None:
    """Run Phase G: Sanity checks."""
    logger.info("\n" + "=" * 70)
    logger.info("PHASE G: SANITY CHECKS")
    logger.info("=" * 70)
    
    results = sanity.check_all_datasets(base_dir, datasets_config)
    sanity.print_check_report(results)


def main():
    parser = argparse.ArgumentParser(
        description="Run MAVE pipeline phases"
    )
    parser.add_argument(
        "--phase",
        choices=["ingest", "normalize", "features", "eval", "all"],
        default="all",
        help="Which phase to run"
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Run sanity checks"
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/mave_datasets.yaml"),
        help="Path to datasets config"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data_processed/mave"),
        help="Base output directory"
    )
    parser.add_argument(
        "-k",
        type=int,
        nargs="+",
        default=[30, 50],
        help="K values for selection (default: 30 50)"
    )
    
    args = parser.parse_args()
    
    logger.info("MAVE Pipeline")
    logger.info(f"Config: {args.config}")
    logger.info(f"Output: {args.output}")
    
    # Ensure directories exist
    args.output.mkdir(parents=True, exist_ok=True)
    
    if args.phase in ["ingest", "all"]:
        run_phase_ingest(args.config, args.output)
    
    if args.phase in ["normalize", "all"]:
        run_phase_normalize(args.config, args.output)
    
    if args.phase in ["features", "all"]:
        run_phase_features(args.output, args.config)
    
    if args.phase in ["eval", "all"]:
        run_phase_evaluation(args.output, args.k)
    
    if args.check:
        run_sanity_checks(args.output, args.config)


if __name__ == "__main__":
    main()

