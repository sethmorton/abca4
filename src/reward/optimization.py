#!/usr/bin/env python3
"""
Optimization routines for ABCA4 variant selection.

Moved from notebooks/03_optimization_dashboard.py to keep notebooks thin.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional
from src.config import logger


class VariantOptimizer:
    """Optimization routines for variant selection with coverage constraints."""

    @staticmethod
    def select_greedy(df: pd.DataFrame, k: int, lambda_penalty: float,
                     cluster_targets: Optional[Dict] = None) -> pd.DataFrame:
        """
        Greedy selection with coverage constraints and bonus for crossing thresholds.

        Args:
            df: DataFrame with variants
            k: Number of variants to select
            lambda_penalty: Coverage penalty weight
            cluster_targets: Dict mapping cluster_id to target coverage threshold

        Returns:
            Selected variants DataFrame
        """
        df_work = df.copy()
        score_col = 'impact_score' if 'impact_score' in df_work.columns else 'model_score'
        if score_col not in df_work.columns:
            raise ValueError("No impact_score or model_score available for optimization")

        # Initialize cluster coverage tracking
        cov = {cid: 0.0 for cid in (cluster_targets.keys() if cluster_targets else [])}

        def _penalty(cov_dict):
            if not cluster_targets:
                return 0.0
            return sum(max(0.0, cluster_targets.get(cid, 0.0) - cov_dict.get(cid, 0.0))
                      for cid in cluster_targets)

        selected = []
        remaining_idx = set(df_work.index)

        for _ in range(k):
            best_idx = None
            best_gain = -np.inf

            for idx in list(remaining_idx):
                row = df_work.loc[idx]
                cid = row.get('cluster_id', None)
                score = float(row.get(score_col, 0.0))

                cov_new = cov.copy()
                coverage_bonus = 0.0

                if cid in (cluster_targets.keys() if cluster_targets else []):
                    old_cov = cov.get(cid, 0.0)
                    new_cov = max(old_cov, score)
                    cov_new[cid] = new_cov
                    # Add small coverage bonus when cluster crosses its threshold (beta=0.05)
                    if new_cov >= cluster_targets[cid] and old_cov < cluster_targets[cid]:
                        coverage_bonus = 0.05

                gain = score - lambda_penalty * (_penalty(cov_new) - _penalty(cov)) + coverage_bonus

                if gain > best_gain:
                    best_gain = gain
                    best_idx = idx

            if best_idx is None:
                break

            remaining_idx.remove(best_idx)
            row = df_work.loc[best_idx]
            cid = row.get('cluster_id', None)

            if cid in (cluster_targets.keys() if cluster_targets else []):
                score_val = float(row.get(score_col, 0.0))
                cov[cid] = max(cov.get(cid, 0.0), score_val)

            selected.append(best_idx)

        # Log coverage shortfall per cluster
        logger.info(f"Optimization complete: λ={lambda_penalty}, K={k}")
        for cid, current_cov in cov.items():
            target = cluster_targets.get(cid, 0.0) if cluster_targets else 0.0
            shortfall = max(0.0, target - current_cov)
            status = "✅ MET" if shortfall == 0 else f"❌ shortfall {shortfall:.3f}"
            logger.info(f"  Cluster {cid}: τⱼ={target:.3f}, cov={current_cov:.3f}, {status}")

        df_sel = df_work.loc[selected].copy()
        df_sel['optimization_score'] = df_sel[score_col]
        df_sel['selected'] = True
        df_sel['rank'] = range(1, len(df_sel) + 1)

        return df_sel

    @staticmethod
    def select_cem(df: pd.DataFrame, k: int, lambda_penalty: float,
                  iters: int = 300, pop: int = 120) -> pd.DataFrame:
        """
        Cross-Entropy Method optimization for variant selection.
        Moved from notebook for reusability.
        """
        score_col = 'impact_score' if 'impact_score' in df.columns else 'model_score'
        if score_col not in df.columns:
            raise ValueError("No impact_score or model_score available for optimization")

        # Precompute cluster targets
        cluster_targets = {}
        if 'cluster_id' in df.columns and 'cluster_target' in df.columns:
            cluster_targets = (
                df[['cluster_id', 'cluster_target']]
                .dropna()
                .drop_duplicates('cluster_id')
                .set_index('cluster_id')['cluster_target']
                .to_dict()
            )

        def _reward(indices, df_subset, score_col, lambda_penalty, cluster_targets):
            sel_df = df_subset.iloc[indices]
            total_score = sel_df[score_col].sum()

            # Coverage penalty
            if cluster_targets:
                cluster_cov = {}
                for idx in indices:
                    row = df_subset.iloc[idx]
                    cid = row.get('cluster_id')
                    if cid in cluster_targets:
                        score = float(row.get(score_col, 0.0))
                        cluster_cov[cid] = max(cluster_cov.get(cid, 0.0), score)

                coverage_penalty = sum(max(0.0, cluster_targets.get(cid, 0.0) - cluster_cov.get(cid, 0.0))
                                     for cid in cluster_targets)
            else:
                coverage_penalty = 0.0

            return total_score - lambda_penalty * coverage_penalty

        n = len(df)
        probs = np.ones(n) / n
        elite_frac = 0.2
        elite_size = max(1, int(pop * elite_frac))

        best_idx = None
        best_reward = -np.inf

        for _ in range(iters):
            samples = []
            rewards = []
            for _ in range(pop):
                idx = np.random.choice(n, size=k, replace=False, p=probs)
                r = _reward(idx, df, score_col, lambda_penalty, cluster_targets)
                samples.append(idx)
                rewards.append(r)

            rewards = np.array(rewards)
            elite_idx = rewards.argsort()[-elite_size:]
            elite_samples = [samples[i] for i in elite_idx]

            flat = np.concatenate(elite_samples)
            counts = np.bincount(flat, minlength=n) + 1e-6  # laplace smoothing
            probs = counts / counts.sum()

            top_elite = elite_samples[-1]
            top_reward = rewards[elite_idx][-1]
            if top_reward > best_reward:
                best_reward = top_reward
                best_idx = top_elite

        df_sel = df.iloc[best_idx].copy()
        df_sel['optimization_score'] = df_sel[score_col]
        df_sel['selected'] = True
        df_sel['rank'] = range(1, len(df_sel) + 1)

        return df_sel

    @staticmethod
    def compute_comparison_metrics(selected_variants: pd.DataFrame,
                                 top_k_variants: pd.DataFrame) -> Dict[str, float]:
        """
        Compute comparison metrics between coverage-aware selection and top-K baseline.

        Returns dict with comparison statistics.
        """
        score_col = 'impact_score' if 'impact_score' in selected_variants.columns else 'model_score'

        selected_total_impact = selected_variants[score_col].sum()
        top_k_total_impact = top_k_variants[score_col].sum()

        impact_diff_pct = ((selected_total_impact - top_k_total_impact) / top_k_total_impact) * 100

        return {
            'selected_total_impact': selected_total_impact,
            'top_k_total_impact': top_k_total_impact,
            'impact_diff_pct': impact_diff_pct,
            'selected_count': len(selected_variants),
            'top_k_count': len(top_k_variants)
        }

    @staticmethod
    def prepare_optimization_data(df_scored: pd.DataFrame,
                                cluster_targets: Optional[Dict] = None) -> pd.DataFrame:
        """
        Prepare DataFrame for optimization by adding cluster_id and cluster_target columns.

        Args:
            df_scored: Scored variants DataFrame
            cluster_targets: Dict mapping cluster names to target thresholds

        Returns:
            DataFrame ready for optimization
        """
        df_opt = df_scored.copy()

        # Add cluster_id column (map cluster names to IDs)
        if 'cluster' in df_opt.columns:
            # Create a simple mapping from cluster names to IDs
            unique_clusters = df_opt['cluster'].unique()
            cluster_id_map = {cluster: f"c{i}" for i, cluster in enumerate(unique_clusters)}
            df_opt['cluster_id'] = df_opt['cluster'].map(cluster_id_map)

            # Add cluster_target column
            if cluster_targets:
                df_opt['cluster_target'] = df_opt['cluster'].map(cluster_targets).fillna(0.0)
            else:
                # Default: use 70th percentile per cluster
                cluster_targets_default = {}
                for cluster_name in unique_clusters:
                    cluster_data = df_opt[df_opt['cluster'] == cluster_name]
                    if not cluster_data.empty:
                        score_col = 'impact_score' if 'impact_score' in cluster_data.columns else 'model_score'
                        if score_col in cluster_data.columns:
                            cluster_targets_default[cluster_name] = min(
                                cluster_data[score_col].quantile(0.7), 0.9
                            )
                        else:
                            cluster_targets_default[cluster_name] = 0.0
                    else:
                        cluster_targets_default[cluster_name] = 0.0

                df_opt['cluster_target'] = df_opt['cluster'].map(cluster_targets_default).fillna(0.0)

        return df_opt
