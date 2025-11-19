#!/usr/bin/env python3
"""Compute regulatory feature set for ABCA4 variants."""

import json
import math
import re
from pathlib import Path
from typing import Dict, Optional

import pandas as pd
import pysam
from pyensembl import EnsemblRelease
from src.config import logger

CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]
CANONICAL_TRANSCRIPT = "ENST00000370225"
ENSEMBL_RELEASE = 109


class RegulatoryFeatureComputer:
    """Derive regulatory heuristics for ABCA4 variants."""

    def __init__(self,
                 annotations_dir: Optional[Path] = None,
                 output_dir: Optional[Path] = None):
        processed_root = CAMPAIGN_ROOT / "data_processed"
        self.annotations_dir = annotations_dir or (processed_root / "annotations")
        self.output_dir = output_dir or (processed_root / "features")
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.gnomad_dir = CAMPAIGN_ROOT / "data_raw" / "gnomad"
        self.domain_config = CAMPAIGN_ROOT / "src" / "data" / "domains" / "abca4_domains.json"

        self.ensembl = EnsemblRelease(ENSEMBL_RELEASE)
        self._ensure_ensembl_ready()
        self.transcript = self.ensembl.transcript_by_id(CANONICAL_TRANSCRIPT)

    def _ensure_ensembl_ready(self) -> None:
        logger.info("Preparing Ensembl cache (%s)...", ENSEMBL_RELEASE)
        try:
            self.ensembl.download()
            self.ensembl.index()
        except Exception as exc:
            logger.warning("Unable to refresh Ensembl cache: %s", exc)

    def load_variants(self) -> Optional[pd.DataFrame]:
        path = self.annotations_dir / "abca4_vus_annotated.parquet"
        if not path.exists():
            logger.error("Annotated variants not found at %s", path)
            return None

        try:
            df = pd.read_parquet(path)
            logger.info("Loaded %s annotated variants", len(df))
            return df
        except Exception as exc:
            logger.error("Unable to read %s: %s", path, exc)
            return None

    def load_domain_config(self) -> Dict:
        with open(self.domain_config, 'r') as handle:
            config = json.load(handle)
        return config

    def _normalize_chrom(self, value: str) -> str:
        return value.replace('chr', '')

    def _variant_key(self, chrom: str, pos: int, ref: str, alt: str) -> str:
        norm_chrom = self._normalize_chrom(chrom)
        return f"{norm_chrom}_{pos}_{ref}_{alt}"

    def _extract_info_value(self, info, field: str, alt_index: int) -> Optional[float]:
        raw = info.get(field)
        if raw is None:
            return None
        if isinstance(raw, (int, float)):
            return float(raw)
        try:
            return float(raw[alt_index])
        except (TypeError, IndexError):
            return None

    def load_gnomad_metrics(self) -> Dict[str, Dict[str, float]]:
        metrics: Dict[str, Dict[str, float]] = {}
        files = [
            (self.gnomad_dir / "gnomad_v4.1_abca4_genome.vcf", "gnomad_genome"),
            (self.gnomad_dir / "gnomad_v4.1_abca4_exome.vcf", "gnomad_exome"),
        ]

        for path, prefix in files:
            if not path.exists():
                logger.warning("Missing gnomAD file: %s", path)
                continue

            logger.info("Parsing %s", path.name)
            try:
                vcf = pysam.VariantFile(path)
            except Exception as exc:
                logger.error("Unable to open %s: %s", path, exc)
                continue

            for record in vcf:
                if not record.alts:
                    continue
                for idx, alt in enumerate(record.alts):
                    key = self._variant_key(record.chrom, record.pos, record.ref, alt)
                    bucket = metrics.setdefault(key, {})

                    af = self._extract_info_value(record.info, 'AF', idx)
                    ac = self._extract_info_value(record.info, 'AC', idx)
                    nhomalt = self._extract_info_value(record.info, 'nhomalt', idx)
                    faf95 = self._extract_info_value(record.info, 'faf95', idx)
                    fafmax_raw = record.info.get('fafmax_faf95_max')
                    grpmax = record.info.get('grpmax')

                    if af is not None:
                        bucket[f"{prefix}_af"] = af
                    if ac is not None:
                        bucket[f"{prefix}_ac"] = ac
                    if nhomalt is not None:
                        bucket[f"{prefix}_nhomalt"] = nhomalt
                    if faf95 is not None:
                        bucket[f"{prefix}_faf95"] = faf95
                    if fafmax_raw is not None:
                        if isinstance(fafmax_raw, (list, tuple)):
                            fafmax_val = max(float(x) for x in fafmax_raw if x is not None)
                        else:
                            fafmax_val = float(fafmax_raw)
                        bucket['faf95_max'] = max(bucket.get('faf95_max', 0.0), fafmax_val)
                    if grpmax is not None:
                        if isinstance(grpmax, (list, tuple)):
                            grp_value = grpmax[0]
                        else:
                            grp_value = grpmax
                        bucket['gnomad_max_population'] = grp_value

        logger.info("Loaded gnomAD metrics for %s variant keys", len(metrics))
        return metrics

    def _map_domains(self, variants: pd.DataFrame, config: Dict) -> pd.DataFrame:
        domain_rows = []
        domains = config.get('domains', [])

        def locate_domain(position: Optional[int]) -> Dict[str, Optional[str]]:
            if position is None:
                return {'domain_label': None, 'domain_class': None, 'domain_distance': math.inf}
            for domain in domains:
                if domain['start'] <= position <= domain['end']:
                    distance = min(position - domain['start'], domain['end'] - position)
                    return {
                        'domain_label': domain['name'],
                        'domain_class': domain['class'],
                        'domain_distance': distance,
                    }
            nearest = min(
                (min(abs(position - d['start']), abs(position - d['end'])) for d in domains),
                default=math.inf,
            )
            return {'domain_label': None, 'domain_class': None, 'domain_distance': nearest}

        variants = variants.copy()
        variants['aa_position'] = variants['protein_change'].apply(self._extract_aa_position)
        domain_info = variants['aa_position'].apply(locate_domain).apply(pd.Series)
        variants = pd.concat([variants, domain_info], axis=1)
        variants['in_domain'] = variants['domain_label'].notna().astype(int)
        return variants

    @staticmethod
    def _extract_aa_position(change: Optional[str]) -> Optional[int]:
        if not isinstance(change, str):
            return None
        match = re.search(r"(\d+)", change)
        if match:
            return int(match.group(1))
        return None

    def add_regulatory_context(self, variants: pd.DataFrame, config: Dict) -> pd.DataFrame:
        regions = config.get('regulatory_regions', [])
        variants = variants.copy()
        variants['regulatory_region'] = 'intergenic'
        variants['regulatory_type'] = 'none'
        variants['regulatory_priority'] = 0.0

        for region in regions:
            mask = (
                variants['chrom'].astype(str).str.replace('chr', '') == region['chrom']
            ) & (
                (variants['pos'] >= region['start']) & (variants['pos'] <= region['end'])
            )
            variants.loc[mask, 'regulatory_region'] = region['name']
            variants.loc[mask, 'regulatory_type'] = region.get('type', 'region')
            variants.loc[mask, 'regulatory_priority'] = region.get('priority', 0.5)

        tss = self.transcript.start
        gene_length = self.transcript.end - self.transcript.start
        variants['distance_to_tss'] = (variants['pos'] - tss).abs()
        variants['relative_position'] = (variants['pos'] - tss) / gene_length
        variants['relative_position'] = variants['relative_position'].clip(0, 1)
        variants['tss_window_score'] = (1 - variants['distance_to_tss'] / 50000).clip(lower=0)
        return variants

    def merge_gnomad(self, variants: pd.DataFrame, metrics: Dict[str, Dict[str, float]]) -> pd.DataFrame:
        if not metrics:
            variants['gnomad_genome_af'] = 0.0
            variants['gnomad_exome_af'] = 0.0
            variants['gnomad_max_af'] = 0.0
            variants['gnomad_combined_nhomalt'] = 0
            variants['gnomad_max_population'] = None
            variants['faf95_max'] = None
            return variants

        lookup = pd.DataFrame.from_dict(metrics, orient='index')
        lookup.index.name = 'variant_key'
        lookup = lookup.reset_index()
        variants = variants.copy()
        variants['variant_key'] = (
            variants['chrom'].astype(str).str.replace('chr', '') + '_' +
            variants['pos'].astype(int).astype(str) + '_' +
            variants['ref'].astype(str) + '_' +
            variants['alt'].astype(str)
        )
        merged = variants.merge(lookup, on='variant_key', how='left')
        merged['gnomad_genome_af'] = merged['gnomad_genome_af'].fillna(0.0)
        merged['gnomad_exome_af'] = merged['gnomad_exome_af'].fillna(0.0)
        merged['gnomad_max_af'] = merged[['gnomad_genome_af', 'gnomad_exome_af']].max(axis=1)
        merged['gnomad_combined_nhomalt'] = (
            merged['gnomad_genome_nhomalt'].fillna(0) + merged['gnomad_exome_nhomalt'].fillna(0)
        )
        return merged

    def compute_regulatory_score(self, variants: pd.DataFrame) -> pd.DataFrame:
        variants = variants.copy()
        variants['faf95_max'] = variants['faf95_max'].fillna(0.0)
        rarity_term = 1 - variants['gnomad_max_af'].clip(upper=0.05)
        variants['regulatory_score'] = (
            variants['regulatory_priority'] * 0.6 +
            variants['tss_window_score'] * 0.2 +
            rarity_term * 0.2
        )
        variants['regulatory_score'] = variants['regulatory_score'].clip(0, 1)
        return variants

    def save(self, df: pd.DataFrame) -> bool:
        output_path = self.output_dir / "regulatory_features.parquet"
        try:
            df.to_parquet(output_path, index=False)
            logger.info("Saved %s regulatory feature rows", len(df))
            return True
        except Exception as exc:
            logger.error("Unable to save regulatory features: %s", exc)
            return False

    def run(self) -> bool:
        variants = self.load_variants()
        if variants is None:
            return False

        config = self.load_domain_config()
        gnomad_metrics = self.load_gnomad_metrics()

        variants = self._map_domains(variants, config)
        variants = self.add_regulatory_context(variants, config)
        variants = self.merge_gnomad(variants, gnomad_metrics)
        variants = self.compute_regulatory_score(variants)

        columns = [
            'variant_id', 'chrom', 'pos', 'ref', 'alt',
            'regulatory_region', 'regulatory_type', 'regulatory_priority',
            'distance_to_tss', 'relative_position', 'tss_window_score',
            'gnomad_genome_af', 'gnomad_exome_af', 'gnomad_max_af',
            'gnomad_combined_nhomalt', 'gnomad_max_population', 'faf95_max',
            'domain_label', 'domain_class', 'domain_distance', 'in_domain',
            'regulatory_score'
        ]
        available = [col for col in columns if col in variants.columns]
        feature_df = variants[available].copy()
        feature_df['distance_to_tss'] = feature_df['distance_to_tss'].fillna(-1)
        feature_df['domain_distance'] = feature_df['domain_distance'].replace(math.inf, -1)

        return self.save(feature_df)


def main() -> None:
    computer = RegulatoryFeatureComputer()
    success = computer.run()
    raise SystemExit(0 if success else 1)


if __name__ == "__main__":
    main()
