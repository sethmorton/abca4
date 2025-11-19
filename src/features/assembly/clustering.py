#!/usr/bin/env python3
"""
Add proper Step 5 clustering to variants_scored.parquet.

Ensures every variant has:
- cluster_id: Domain or consequence-based cluster assignment
- cluster_target: τⱼ coverage threshold per cluster
- coverage_by_cluster: cov_j(S) = max impact_score in cluster
"""

import sys
from pathlib import Path
from typing import Optional, Dict
import json

import pandas as pd
from src.config import logger
from src.config import load_gene_config

CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]


class ClusteringProcessor:
    """Add Step 5 clustering to scored variants using config-driven parameters."""

    def __init__(self, gene_name: str, features_dir: Optional[Path] = None, config: Optional[Dict] = None):
        """
        Initialize clustering processor.
        
        Args:
            gene_name: Gene symbol (e.g., "ABCA4") - REQUIRED, no default
            features_dir: Directory for features
            config: Optional pre-loaded gene config. If not provided, loads from config file.
            
        Raises:
            ValueError: If gene_name is empty or invalid
        """
        if not gene_name or not isinstance(gene_name, str) or gene_name.strip() == "":
            raise ValueError("gene_name is required and must be a non-empty string")
        
        self.gene_name = gene_name
        self.config = config or load_gene_config(gene_name)
        
        processed_root = CAMPAIGN_ROOT / "data_processed"
        self.features_dir = features_dir or (processed_root / "features")
        self.features_dir.mkdir(parents=True, exist_ok=True)

    def load_scored_variants(self) -> Optional[pd.DataFrame]:
        """Load variants with scored features (with impact scores and metadata)."""
        # Use gene prefix for all filenames
        file_prefix = self.config.get("file_prefix", self.gene_name.lower())
        
        # First try variants_features_raw.parquet (has impact scores and metadata)
        features_raw_path = self.features_dir / f"{file_prefix}_variants_features_raw.parquet"
        scored_path = self.features_dir / f"{file_prefix}_variants_scored.parquet"
        metadata_path = self.features_dir / f"{file_prefix}_variants_scored.metadata.json"

        load_path = None
        if features_raw_path.exists():
            load_path = features_raw_path
        elif scored_path.exists():
            load_path = scored_path
        else:
            logger.error(f"Missing both variants_features_raw.parquet and variants_scored.parquet")
            return None

        try:
            df = pd.read_parquet(load_path)
            logger.info(f"Loaded {len(df)} scored variants from {load_path.name}")

            # Attempt to restore attrs from sidecar metadata if present
            try:
                if metadata_path.exists():
                    with open(metadata_path, "r", encoding="utf-8") as f:
                        restored_attrs = json.load(f)
                    df.attrs.update(restored_attrs)
                    logger.info("Restored DataFrame attrs from metadata sidecar")
                elif load_path == features_raw_path:
                    # If loading from features_raw and no sidecar, attrs should already be there
                    logger.info(f"Loaded from {load_path.name} with existing attrs: {list(df.attrs.keys())}")
            except Exception as meta_err:
                logger.warning(f"Could not restore attrs metadata: {meta_err}")
            return df
        except Exception as e:
            logger.error(f"Failed to load scored variants: {e}")
            return None

    def apply_domain_boost(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply domain-aware scoring boost to model_score using config factors."""
        df_boosted = df.copy()

        if 'domain' not in df_boosted.columns or 'model_score' not in df_boosted.columns:
            logger.warning("Domain or model_score column missing; skipping domain boost")
            return df_boosted

        # Get domain boost factors from config, with defaults for common non-coding domains
        domain_boost = self.config.get('domain_boost_factors', {})
        
        # Add sensible defaults for non-coding domains if not in config
        defaults = {
            'CTD': 1.0,         # No boost
            'intronic': 1.0,    # No boost
            'utr': 0.95,        # Slight penalty
            'other': 1.0,       # No boost
            'unknown': 1.0,     # No boost
            'regulatory': 1.0,  # No boost
            'synonymous': 0.9,  # Slight penalty
            'splice_region': 1.0, # No boost
        }
        
        for domain, factor in defaults.items():
            if domain not in domain_boost:
                domain_boost[domain] = factor

        original_scores = df_boosted['model_score'].copy()
        df_boosted['model_score'] = df_boosted.apply(
            lambda row: row['model_score'] * domain_boost.get(row['domain'], 1.0),
            axis=1
        ).clip(0, 1)

        boosted_count = (df_boosted['model_score'] != original_scores).sum()
        logger.info(f"Applied domain-aware scoring boost to {boosted_count} variants")
        logger.info(f"Score range: [{original_scores.min():.4f}, {original_scores.max():.4f}] → "
                   f"[{df_boosted['model_score'].min():.4f}, {df_boosted['model_score'].max():.4f}]")

        return df_boosted

    def apply_hotspot_overrides(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply hotspot cluster overrides for known mechanisms."""
        df_overridden = df.copy()

        if 'pos' not in df_overridden.columns or 'domain' not in df_overridden.columns:
            logger.warning("Position or domain column missing; skipping hotspot overrides")
            return df_overridden

        # Hotspot overrides for known Stargardt mechanisms
        hotspot_overrides = {
            # Format: (start_pos, end_pos): override_domain
            # TODO: Fill in actual hotspot ranges based on literature
            # (588, 588): 'NBD1_hotspot',
            # (768, 768): 'TMD1_hotspot',
        }

        if hotspot_overrides:
            for (start_pos, end_pos), override_domain in hotspot_overrides.items():
                pos_mask = (df_overridden['pos'] >= start_pos) & (df_overridden['pos'] <= end_pos)
                df_overridden.loc[pos_mask, 'domain'] = override_domain
                logger.info(f"Applied hotspot override: positions {start_pos}-{end_pos} → {override_domain}")

        return df_overridden

    def assign_clusters_domain(self, df: pd.DataFrame) -> pd.DataFrame:
        """Assign clusters based on domain (primary) with consequence sub-clustering."""
        df_clustered = df.copy()

        if "domain" in df_clustered.columns and df_clustered["domain"].notna().any():
            # Group by domain first
            df_clustered["cluster_id"] = df_clustered["domain"].fillna("unknown")

            # For domains with many variants, sub-cluster by consequence
            large_threshold = self.config.get('clustering', {}).get('large_cluster_threshold', 10)
            domain_counts = df_clustered["cluster_id"].value_counts()
            large_domains = domain_counts[domain_counts > large_threshold].index

            for domain in large_domains:
                mask = df_clustered["cluster_id"] == domain
                if "vep_consequence" in df_clustered.columns:
                    # Create sub-clusters within large domains
                    sub_cluster = df_clustered.loc[mask, "vep_consequence"].fillna("other")
                    df_clustered.loc[mask, "cluster_id"] = sub_cluster.astype(str).radd(f"{domain}_")

            logger.info(f"Domain-based clustering: {df_clustered['cluster_id'].nunique()} clusters")
        else:
            # Fallback to consequence clustering
            logger.warning("Domain info unavailable; falling back to consequence clustering")
            df_clustered = self.assign_clusters_consequence(df_clustered)

        return df_clustered

    def assign_clusters_consequence(self, df: pd.DataFrame) -> pd.DataFrame:
        """Assign clusters based on variant consequence (primary mechanism)."""
        df_clustered = df.copy()

        def consequence_to_cluster(consequence_str):
            """Map VEP consequence to cluster category."""
            if not isinstance(consequence_str, str) or pd.isna(consequence_str):
                return "other"

            consequence_lower = consequence_str.lower()

            # Check for LoF (high priority)
            if any(x in consequence_lower for x in ["frameshift", "stop_gained", "stop_lost"]):
                if "frameshift" in consequence_lower:
                    return "LoF_frameshift"
                return "LoF_stop"

            # Check for canonical splice variants
            if any(x in consequence_lower for x in ["splice_acceptor", "splice_donor"]):
                if "acceptor" in consequence_lower:
                    return "LoF_splice_acceptor"
                return "LoF_splice_donor"

            # Check for missense
            if "missense" in consequence_lower:
                return "missense"

            # Check for splice region/proxy (moderate impact)
            if any(x in consequence_lower for x in ["splice_region", "splice_polypyrimidine"]):
                return "splice_region"

            # Check for indels
            if any(x in consequence_lower for x in ["inframe_deletion", "inframe_insertion", "disruptive"]):
                return "inframe_indel"

            # Check for synonymous
            if "synonymous" in consequence_lower:
                return "synonymous"

            # Check for intron
            if "intron" in consequence_lower:
                return "intron"

            # Check for UTR
            if "utr" in consequence_lower:
                return "utr"

            return "other"

        if "vep_consequence" in df_clustered.columns:
            df_clustered["cluster_id"] = df_clustered["vep_consequence"].apply(consequence_to_cluster)
            logger.info("Using vep_consequence for mechanism-based clustering")
        else:
            df_clustered["cluster_id"] = "unknown"
            logger.warning("vep_consequence not found; using 'unknown' for all")

        return df_clustered

    def compute_cluster_targets_biology_aware(self, df: pd.DataFrame,
                                            clinical_significance_col: str = "clinical_significance") -> pd.DataFrame:
        """Compute biology-aware τⱼ targets: median P/LP score if available, else 70th percentile."""
        df_with_targets = df.copy()

        cluster_targets = {}

        for cluster_name, group in df_with_targets.groupby("cluster_id"):
            # Check for pathogenic/likely pathogenic labels in this cluster
            if clinical_significance_col in group.columns:
                def _is_pathogenic(sig):
                    return "pathogenic" in str(sig).lower() and "benign" not in str(sig).lower()
                pathogenic_mask = group[clinical_significance_col].apply(_is_pathogenic)
                pathogenic_scores = group.loc[pathogenic_mask, "model_score"]

                if not pathogenic_scores.empty:
                    # Use median of P/LP variants
                    tau_j = pathogenic_scores.median()
                    logger.info(f"  {cluster_name}: P/LP median τⱼ={tau_j:.4f}")
                else:
                    # Use 70th percentile of all scores in cluster, cap at 0.9
                    tau_j = min(group["model_score"].quantile(0.7), 0.9)
                    logger.info(f"  {cluster_name}: 70th percentile τⱼ={tau_j:.4f}")
            else:
                # Fallback: 70th percentile
                tau_j = min(group["model_score"].quantile(0.7), 0.9)
                logger.info(f"  {cluster_name}: 70th percentile τⱼ={tau_j:.4f} (no clinical data)")

            cluster_targets[cluster_name] = tau_j

        df_with_targets["cluster_target"] = df_with_targets["cluster_id"].map(cluster_targets)

        logger.info(f"Computed biology-aware targets for {len(cluster_targets)} clusters")
        return df_with_targets

    def compute_cluster_targets(self, df: pd.DataFrame, threshold_factor: float = 0.8) -> pd.DataFrame:
        """Legacy method: Compute τⱼ = threshold_factor × max(model_score) per cluster."""
        df_with_targets = df.copy()

        cluster_targets = {}
        for cluster_id, group in df_with_targets.groupby("cluster_id"):
            max_score = group["model_score"].max()
            tau_j = max_score * threshold_factor
            cluster_targets[cluster_id] = tau_j
            logger.info(f"  {cluster_id}: max_score={max_score:.4f}, τⱼ={tau_j:.4f}")

        df_with_targets["cluster_target"] = df_with_targets["cluster_id"].map(cluster_targets)

        logger.info(f"Computed legacy targets for {len(cluster_targets)} clusters (threshold={threshold_factor})")
        return df_with_targets

    def compute_cluster_coverage(self, df: pd.DataFrame) -> pd.DataFrame:
        """Compute cov_j(S) = max(model_score) in each cluster."""
        df_with_coverage = df.copy()
        
        cluster_coverage = {}
        for cluster_id, group in df_with_coverage.groupby("cluster_id"):
            max_score = group["model_score"].max()
            cluster_coverage[cluster_id] = max_score
            logger.info(f"  {cluster_id}: cov_j(S)={max_score:.4f}")
        
        df_with_coverage["coverage_by_cluster"] = df_with_coverage["cluster_id"].map(cluster_coverage)
        
        logger.info(f"Computed coverage metrics for {len(cluster_coverage)} clusters")
        return df_with_coverage

    def save_clustered_variants(self, df: pd.DataFrame) -> bool:
        """Save variants with clustering info."""
        file_prefix = self.config.get("file_prefix", self.gene_name.lower())
        output_path = self.features_dir / f"{file_prefix}_variants_scored.parquet"
        metadata_path = self.features_dir / f"{file_prefix}_variants_scored.metadata.json"
        
        try:
            df.to_parquet(output_path, index=False)
            logger.info(f"Saved {len(df)} clustered variants to {output_path}")
            # Persist DataFrame.attrs to sidecar metadata for reproducibility/traceability
            try:
                # Ensure attrs are JSON-serializable; fallback to string for non-serializable items
                serializable_attrs = {}
                for k, v in df.attrs.items():
                    try:
                        json.dumps(v)
                        serializable_attrs[k] = v
                    except Exception:
                        serializable_attrs[k] = str(v)
                with open(metadata_path, "w", encoding="utf-8") as f:
                    json.dump(serializable_attrs, f, indent=2)
                logger.info(f"Wrote DataFrame attrs to {metadata_path}")
            except Exception as meta_err:
                logger.warning(f"Failed to write attrs metadata: {meta_err}")
            
            # Log summary
            if "cluster_id" in df.columns:
                cluster_counts = df["cluster_id"].value_counts()
                logger.info(f"Cluster distribution:\n{cluster_counts}")
            
            return True
        except Exception as e:
            logger.error(f"Failed to save clustered variants: {e}")
            return False

    def run(self, clustering_mode: str = "domain", use_biology_aware_targets: bool = True) -> bool:
        """Execute full clustering pipeline with domain-aware processing."""
        logger.info("Starting clustering processor...")

        # Load
        df = self.load_scored_variants()
        if df is None:
            return False
        # Keep attrs snapshot to propagate through transformations
        _attrs = dict(df.attrs)

        # STRICT VALIDATION: Required columns must be present
        required_columns = ["model_score", "variant_id"]
        if clustering_mode == "domain":
            required_columns.append("domain")
        else:
            required_columns.append("vep_consequence")

        missing_columns = [col for col in required_columns if col not in df.columns]

        if missing_columns:
            error_msg = f"""
❌ CLUSTERING FAILED: Missing required columns.

Required columns: {required_columns}
Missing columns: {missing_columns}
Available columns: {list(df.columns)}

ACTION REQUIRED:
Ensure variants_scored.parquet contains required columns.
Run feature engineering pipeline first if these are missing.
            """.strip()
            logger.error(error_msg)
            raise SystemExit(1)

        logger.info(f"✅ Required columns present: {required_columns}")

        # Apply domain boost if domain column available
        if "domain" in df.columns:
            logger.info("Applying domain-aware scoring boost...")
            new_df = self.apply_domain_boost(df)
            new_df.attrs.update(_attrs)
            df = new_df

        # Apply hotspot overrides
        logger.info("Applying hotspot overrides...")
        new_df = self.apply_hotspot_overrides(df)
        new_df.attrs.update(df.attrs if df is not None else _attrs)
        df = new_df

        # Assign clusters
        logger.info(f"Assigning clusters using {clustering_mode} mode...")
        if clustering_mode == "domain":
            new_df = self.assign_clusters_domain(df)
        else:
            new_df = self.assign_clusters_consequence(df)
        new_df.attrs.update(df.attrs)
        df = new_df

        # Compute targets (biology-aware or legacy)
        logger.info("Computing cluster targets...")
        if use_biology_aware_targets:
            new_df = self.compute_cluster_targets_biology_aware(df)
        else:
            new_df = self.compute_cluster_targets(df, threshold_factor=0.8)
        new_df.attrs.update(df.attrs)
        df = new_df

        # Compute coverage
        logger.info("Computing cluster coverage...")
        new_df = self.compute_cluster_coverage(df)
        new_df.attrs.update(df.attrs)
        df = new_df

        # Save
        logger.info("Saving results...")
        return self.save_clustered_variants(df)


def main():
    """Main entry point - requires gene to be specified via --gene or GENE_NAME env var."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Add clustering to scored variants")
    parser.add_argument("--gene", type=str, default=os.getenv("GENE_NAME", None),
                       help="Gene symbol (required: pass via --gene or GENE_NAME env var)")
    args = parser.parse_args()
    
    if not args.gene:
        parser.error(
            "❌ ERROR: Gene symbol required but not provided.\n"
            "Please specify one of:\n"
            "  1. Command line: python script.py --gene GENE_NAME\n"
            "  2. Environment: export GENE_NAME=GENE_NAME"
        )
    
    try:
        processor = ClusteringProcessor(args.gene)
        success = processor.run()
        sys.exit(0 if success else 1)
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    import os
    main()
