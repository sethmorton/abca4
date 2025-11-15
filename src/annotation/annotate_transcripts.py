#!/usr/bin/env python3
"""
Add transcript and functional annotations to ABCA4 variants.

Uses VEP or pyensembl to add transcript IDs, protein changes, exon distances, etc.
"""

import pandas as pd
from pathlib import Path
import sys
from typing import Optional, Dict, Any, List
import requests
from pyensembl import EnsemblRelease
from src.config import logger

CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]

CANONICAL_TRANSCRIPT = "ENST00000370225"
ENSEMBL_RELEASE = 109  # GRCh38

class VariantAnnotator:
    """Add transcript and functional annotations to variants."""

    def __init__(self, input_dir: Optional[Path] = None, output_dir: Optional[Path] = None,
                 ensembl_release: int = ENSEMBL_RELEASE):
        data_root = CAMPAIGN_ROOT / "data_processed"
        self.input_dir = input_dir or (data_root / "variants")
        self.output_dir = output_dir or (data_root / "annotations")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.ensembl = EnsemblRelease(ensembl_release)
        self._ensure_ensembl_ready()
        self.transcript = self.ensembl.transcript_by_id(CANONICAL_TRANSCRIPT)

    def fetch_vep_annotations(self, variants_df: pd.DataFrame) -> Dict[str, Dict[str, Any]]:
        """Query Ensembl VEP REST API for HGVS + consequence annotations."""
        url = "https://rest.ensembl.org/vep/human/region"
        headers = {"Content-Type": "application/json", "Accept": "application/json"}
        annotations: Dict[str, Dict[str, Any]] = {}
        batch_size = 100

        rows = variants_df[['variant_id', 'chrom', 'pos', 'ref', 'alt']].to_dict(orient='records')
        for i in range(0, len(rows), batch_size):
            batch = rows[i:i + batch_size]
            payload_variants = []
            lookup: Dict[str, str] = {}
            for row in batch:
                start = int(row['pos'])
                end = start
                allele = f"{row['ref']}/{row['alt']}"
                variant_str = f"{row['chrom']} {start} {end} {allele}"
                payload_variants.append(variant_str)
                lookup[variant_str] = row['variant_id']

            resp = requests.post(url, headers=headers,
                                 json={"variants": payload_variants, "canonical": 1, "protein": 1, "hgvs": 1},
                                 timeout=60)
            if resp.status_code != 200:
                logger.warning("VEP request failed with status %s: %s", resp.status_code, resp.text[:200])
                continue

            for item in resp.json():
                variant_str = item.get('input')
                variant_id = lookup.get(variant_str)
                if not variant_id:
                    continue

                transcript_data = None
                for consequence in item.get('transcript_consequences', []):
                    if consequence.get('transcript_id') == CANONICAL_TRANSCRIPT:
                        transcript_data = consequence
                        break

                if not transcript_data:
                    continue

                annotations[variant_id] = {
                    'transcript_id': transcript_data.get('transcript_id'),
                    'hgvs_c': transcript_data.get('hgvsc'),
                    'hgvs_p': transcript_data.get('hgvsp'),
                    'protein_change': transcript_data.get('hgvsp'),
                    'consequence_terms': ','.join(transcript_data.get('consequence_terms', [])),
                    'impact': transcript_data.get('impact'),
                    'exon_number': transcript_data.get('exon'),
                    'intron_number': transcript_data.get('intron'),
                }

        return annotations

    def _ensure_ensembl_ready(self) -> None:
        logger.info("Ensuring Ensembl release %s metadata is available...", self.ensembl.release)
        try:
            self.ensembl.download()
            self.ensembl.index()
        except Exception as exc:
            logger.warning("Ensembl setup encountered an issue: %s", exc)

    def load_variants(self) -> Optional[pd.DataFrame]:
        """Load filtered ABCA4 variants."""
        variants_path = self.input_dir / "abca4_clinvar_vus.parquet"

        if not variants_path.exists():
            logger.error(f"Variants file not found: {variants_path}")
            return None

        logger.info(f"Loading variants from {variants_path}")
        try:
            df = pd.read_parquet(variants_path)
            logger.info(f"Loaded {len(df)} variants")
            return df
        except Exception as e:
            logger.error(f"Failed to load variants: {e}")
            return None

    def add_transcript_annotations(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Add transcript-level annotations via VEP + pyensembl."""
        logger.info("Adding transcript annotations via Ensembl APIs...")

        annotations = self.fetch_vep_annotations(variants_df)
        annotated_rows: List[Dict[str, Any]] = []

        for _, row in variants_df.iterrows():
            record = row.to_dict()
            vep = annotations.get(record['variant_id'], {})
            record['transcript_id'] = vep.get('transcript_id', CANONICAL_TRANSCRIPT)
            record['protein_change'] = vep.get('protein_change')
            record['hgvs_c'] = vep.get('hgvs_c')
            record['hgvs_p'] = vep.get('hgvs_p')
            record['vep_consequence'] = vep.get('consequence_terms')
            record['vep_impact'] = vep.get('impact')

            exon_ctx = self._compute_structural_context(int(record['pos']))
            record.update(exon_ctx)

            annotated_rows.append(record)

        annotated_df = pd.DataFrame(annotated_rows)
        logger.info("Annotated %s variants with VEP and pyensembl context", len(annotated_df))
        return annotated_df

    def _compute_structural_context(self, position: int) -> Dict[str, Any]:
        region = 'intergenic'
        exon_number = None
        intron_distance = None

        if self.transcript.start <= position <= self.transcript.end:
            region = 'intronic'
            nearest_distance = None

            for exon_number, exon in enumerate(self.transcript.exons, start=1):
                if exon.start <= position <= exon.end:
                    region = 'exonic'
                    intron_distance = min(position - exon.start, exon.end - position)
                    break
                else:
                    distance = min(abs(position - exon.start), abs(position - exon.end))
                    if nearest_distance is None or distance < nearest_distance:
                        nearest_distance = distance

            if region != 'exonic' and nearest_distance is not None:
                intron_distance = nearest_distance
        else:
            if position < self.transcript.start:
                region = 'upstream'
                intron_distance = self.transcript.start - position
            else:
                region = 'downstream'
                intron_distance = position - self.transcript.end

        coding_class = 'noncoding'
        if region == 'exonic':
            coding_class = 'coding_snv'
        elif region in {'intronic', 'upstream', 'downstream'}:
            coding_class = region

        return {
            'genomic_region': region,
            'exon_number': exon_number,
            'intron_distance': intron_distance,
            'coding_impact': coding_class
        }

    def add_genomic_annotations(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Add genomic context annotations."""
        logger.info("Adding genomic context annotations...")

        annotated_df = variants_df.copy()
        annotated_df['gene_start'] = self.transcript.start
        annotated_df['gene_end'] = self.transcript.end
        annotated_df['distance_to_gene_start'] = annotated_df['pos'] - self.transcript.start
        annotated_df['distance_to_gene_end'] = self.transcript.end - annotated_df['pos']
        annotated_df['within_gene_span'] = (
            (annotated_df['pos'] >= self.transcript.start) &
            (annotated_df['pos'] <= self.transcript.end)
        )

        if 'genomic_region' not in annotated_df.columns:
            annotated_df['genomic_region'] = annotated_df['within_gene_span'].map(
                {True: 'intronic', False: 'intergenic'}
            )

        logger.info("Added genomic context annotations")
        return annotated_df

    def validate_annotations(self, annotated_df: pd.DataFrame) -> bool:
        """Validate annotation quality."""
        logger.info("Validating annotations...")

        # Check for required columns
        required_cols = ['transcript_id', 'genomic_region', 'coding_impact']
        missing_cols = [col for col in required_cols if col not in annotated_df.columns]

        if missing_cols:
            logger.error(f"Missing required annotation columns: {missing_cols}")
            return False

        # Check annotation completeness
        completeness = annotated_df[required_cols].notna().mean()
        logger.info(f"Annotation completeness: {completeness}")

        if completeness.min() < 0.8:
            logger.warning("Low annotation completeness detected")

        return True

    def save_annotated_variants(self, annotated_df: pd.DataFrame) -> bool:
        """Save annotated variants."""
        output_path = self.output_dir / "abca4_vus_annotated.parquet"

        try:
            annotated_df.to_parquet(output_path, index=False)
            logger.info(f"Saved {len(annotated_df)} annotated variants to {output_path}")

            # Also save summary CSV
            csv_path = self.output_dir / "abca4_vus_annotated.csv"
            annotated_df.to_csv(csv_path, index=False)
            logger.info(f"Also saved summary CSV: {csv_path}")

            return True
        except Exception as e:
            logger.error(f"Failed to save annotated variants: {e}")
            return False

    def run(self) -> bool:
        """Run the complete annotation process."""
        logger.info("Starting variant annotation process...")

        # Load variants
        variants_df = self.load_variants()
        if variants_df is None:
            return False

        # Add transcript annotations
        annotated_df = self.add_transcript_annotations(variants_df)

        # Add genomic annotations
        annotated_df = self.add_genomic_annotations(annotated_df)

        # Validate annotations
        if not self.validate_annotations(annotated_df):
            logger.error("Annotation validation failed")
            return False

        # Save results
        if not self.save_annotated_variants(annotated_df):
            return False

        logger.info("Variant annotation process completed successfully!")
        return True

def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(description="Annotate ABCA4 variants with transcript context")
    parser.add_argument("--input-dir", type=Path,
                        help="Input variants directory (default: campaigns/abca4/data_processed/variants)")
    parser.add_argument("--output-dir", type=Path,
                        help="Output annotations directory (default: campaigns/abca4/data_processed/annotations)")
    parser.add_argument("--ensembl-release", type=int, default=ENSEMBL_RELEASE,
                        help="Ensembl release to use (default: 109)")

    args = parser.parse_args()

    annotator = VariantAnnotator(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        ensembl_release=args.ensembl_release,
    )
    success = annotator.run()
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
