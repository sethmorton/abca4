#!/usr/bin/env python3
"""
Documentation strings and markdown templates for ABCA4 pipeline notebooks.

Moved from notebooks to keep them thin and consolidate duplicate content.
"""

from typing import Dict, Any


class PipelineDocs:
    """Documentation strings and templates for the ABCA4 pipeline."""

    @staticmethod
    def get_feature_engineering_intro() -> str:
        """Introduction markdown for feature engineering notebook."""
        return """
# ABCA4 Variant Analysis: Feature Engineering

This notebook loads and processes genetic variant features for ABCA4, the gene associated with Stargardt disease.

## Pipeline Overview

1. **Load annotated variants** from ClinVar VUS dataset
2. **Compute/combine features** from multiple sources:
   - AlphaMissense (missense prediction)
   - SpliceAI (splice disruption)
   - Conservation scores (phyloP)
   - Population frequencies (gnomAD)
   - Protein domain assignments
3. **Quality validation** and null rate checking
4. **Impact score computation** with domain-aware weighting
5. **Export** for downstream optimization

## Key Outputs

- `variants_features_raw.parquet`: All features with impact scores
- `variants_scored.parquet`: Clustered variants with coverage targets
- QC plots and validation metrics
"""

    @staticmethod
    def get_data_quality_audit() -> str:
        """Data quality audit section."""
        return """
#### Feature Quality Audit

Validating data integrity and coverage across all feature sources.
"""

    @staticmethod
    def get_impact_scoring_intro() -> str:
        """Impact scoring section introduction."""
        return """
## Step 4: Impact Score Construction

Combine multiple signals into a unified impact score for variant prioritization.

### Scoring Formula

```
impact_score = 0.6 × model_score +
               0.2 × cons_scaled +
               0.1 × domain_flag +
               0.1 × splice_prox_flag -
               0.3 × af_penalty
```

### Component Explanations

- **model_score**: AlphaMissense/SpliceAI/conservation weighted combination
- **cons_scaled**: Normalized phyloP conservation score (0-1)
- **domain_flag**: 1 if variant falls in defined protein domain
- **splice_prox_flag**: 1 if splice-related consequence or near exon boundary
- **af_penalty**: Monotonic penalty for common variants (rarity bonus)
"""

    @staticmethod
    def get_clustering_intro() -> str:
        """Clustering section introduction."""
        return """
## Step 5: Clustering & Coverage Targets

Group variants by biological mechanism and compute coverage thresholds.

### Clustering Strategy

- **Primary grouping**: Protein domain (NBD1, NBD2, TMD1, TMD2, CTD)
- **Sub-clustering**: By consequence within large domains
- **Fallback**: Consequence-based clustering if domains unavailable

### Coverage Targets (τⱼ)

Biology-aware thresholds computed as:
- **With clinical data**: Median impact_score of P/LP variants in cluster
- **Without clinical data**: 70th percentile of cluster scores (capped at 0.9)
"""

    @staticmethod
    def get_optimization_intro() -> str:
        """Optimization dashboard introduction."""
        return """
# ABCA4 Variant Selection: Optimization Dashboard

Interactive optimization of variant panels with coverage constraints.

## Selection Strategies

- **Greedy Coverage**: Maximizes total impact while meeting cluster thresholds
- **Top-K by Score**: Simple ranking by impact score (baseline)
- **CEM**: Cross-Entropy Method for complex optimization
- **Random**: Random sampling (control)

## Coverage Constraints

Each cluster j has target threshold τⱼ. Selection optimizes:

```
max Σᵢ scoreᵢ - λ × Σⱼ max(0, τⱼ - covⱼ(S))
```

Where covⱼ(S) = max score of selected variants in cluster j.
"""

    @staticmethod
    def get_comparison_panel_intro() -> str:
        """Comparison panel documentation."""
        return """
#### Comparison Panel: Coverage-aware vs Top-K by Score

Side-by-side evaluation of optimization strategies:

- **Coverage-aware**: Greedy selection with cluster constraints
- **Top-K by Score**: Baseline ranking by impact score only

**Success Criteria:**
- Total impact within 5% of top-K baseline
- Improved cluster coverage distribution
- No catastrophic impact loss
"""

    @staticmethod
    def get_cro_planning_intro() -> str:
        """CRO planning introduction."""
        return """
## Experimental Planning & CRO Coordination

Translate selected variants into experimental work packages.

### Assay Mapping

- **DSF/SEC**: Protein stability and oligomerization
- **Functional**: Transport activity assays
- **Trafficking**: Protein localization studies
- **RNA-seq**: Expression and splicing analysis
- **Splice Minigene**: In vitro splicing assays

### Work Package Generation

- Automatic assay assignment based on variant mechanisms
- Sample count optimization
- Timeline and resource estimation
"""

    @staticmethod
    def get_feature_status_template(feature_sources: Dict[str, Any]) -> str:
        """Template for feature loading status display."""
        status_lines = []
        for name, info in feature_sources.items():
            status_lines.append(f"  ✅ {name}: {info.get('rows', 0)} rows, {info.get('columns', 0)} columns")

        return f"""
📊 Feature Loading Status:
{chr(10).join(status_lines)}

📄 Status written to: data_processed/features/features_status.json
"""

    @staticmethod
    def get_merge_validation_template(merge_validation: list) -> str:
        """Template for merge validation display."""
        if not merge_validation:
            return "No merge validation data available."

        rows = []
        for item in merge_validation:
            rows.append([
                item.get('feature_source', 'N/A'),
                str(item.get('duplicates_removed', 0)),
                str(item.get('rows_before_merge', 0)),
                str(item.get('rows_after_merge', 0))
            ])

        # Simple table format
        header = "| Feature Source | Duplicates Removed | Rows Before | Rows After |"
        separator = "|---------------|-------------------|-------------|------------|"
        table_rows = [f"| {' | '.join(row)} |" for row in rows]

        return "\n".join([header, separator] + table_rows)

    @staticmethod
    def get_null_rates_template(null_rates: Dict[str, Any]) -> str:
        """Template for null rates display."""
        if not null_rates:
            return "No null rates data available."

        rows = []
        for col, data in null_rates.items():
            rows.append([
                col,
                data.get('description', 'N/A'),
                str(data.get('non_null_count', 0)),
                f"{data.get('non_null_pct', 0):.1f}%"
            ])

        # Simple table format
        header = "| Column | Feature | Non-Null Count | Non-Null % |"
        separator = "|--------|---------|----------------|-------------|"
        table_rows = [f"| {' | '.join(row)} |" for row in rows]

        return "\n".join([header, separator] + table_rows)

    @staticmethod
    def get_cluster_summary_template(df_clusters: Any) -> str:
        """Template for cluster summary display."""
        try:
            total_clusters = df_clusters['cluster_id'].nunique()
            cluster_counts = df_clusters['cluster_id'].value_counts()
            top_clusters = cluster_counts.head().to_dict()

            summary = f"""
**Total clusters:** {total_clusters}

**Top clusters by size:**
"""

            for cluster, count in top_clusters.items():
                summary += f"- {cluster}: {count} variants\n"

            return summary
        except Exception:
            return "Cluster summary not available."

    @staticmethod
    def get_optimization_results_template(results: Dict[str, Any]) -> str:
        """Template for optimization results display."""
        strategy = results.get('strategy', 'Unknown')
        k = results.get('k', 0)
        selected_count = len(results.get('selected_variants', []))

        template = f"""
### Results

- Strategy: {strategy}
- K: {k}
- Selected: {selected_count}
"""

        # Add comparison metrics if available
        top_k = results.get('top_k_by_score')
        selected = results.get('selected_variants')

        if top_k is not None and selected is not None:
            score_col = 'impact_score' if 'impact_score' in selected.columns else 'model_score'
            selected_total = selected[score_col].sum()
            top_k_total = top_k[score_col].sum()
            impact_diff_pct = ((selected_total - top_k_total) / top_k_total) * 100

            template += f"""

#### Comparison Panel: Coverage-aware vs Top-K by Score

| Metric | Coverage-aware | Top-K by Score | Difference |
|--------|---------------|----------------|------------|
| Total Impact | {selected_total:.3f} | {top_k_total:.3f} | {impact_diff_pct:+.1f}% |
| Variants | {len(selected)} | {len(top_k)} | - |
"""

            if impact_diff_pct < -5:
                template += "\n⚠️ **Warning:** Coverage-aware result has >5% lower total impact than top-K baseline."
            else:
                template += "\n✅ Coverage-aware result maintains similar total impact to top-K baseline."

        return template
