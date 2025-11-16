# generate_protocols.py
import os
from pathlib import Path
import pandas as pd
import numpy as np
from datetime import datetime
from groq import Groq               # <-- NEW
from config import get_key         # <-- same as before

# === CONFIG ===

NOTEBOOKS_DIR = Path(__file__).resolve().parent
INPUT_PARQUET = NOTEBOOKS_DIR / ".datasets" / "variants_scored1.parquet"

OUTPUT_DIR = NOTEBOOKS_DIR / "protocols"
TOP_N = 8                         # Number of variants to test (budget)
MODEL = "llama-3.3-70b-versatile" # <-- FREE Groq model (or "gemma2-9b-it")
TEMPERATURE = 0.3

# Groq client (free tier)
client = Groq(api_key=get_key())

# Create output dirs
Path(OUTPUT_DIR).mkdir(exist_ok=True)

print(f"Loading variants from {INPUT_PARQUET}...")
if not os.path.exists(INPUT_PARQUET):
    raise FileNotFoundError(f"Parquet file not found: {INPUT_PARQUET}")

# === 1. Load Parquet ===
df = pd.read_parquet(INPUT_PARQUET)
print(f"Loaded {len(df)} variants. Columns: {list(df.columns)}")

# === 2. Compute Blended Impact Score ===
SCORE_COLUMNS = {
    'alphamissense_score': 0.5,
    'spliceai_max_score': 0.2,
    'phylop_score': 0.2,
    'phastcons_score': 0.2,
    'conservation_score': 0.2,
    'tss_window_score': 0.2,
    'regulatory_score': 0.2,
    'model_score': 0.2
}

available = {col: w for col, w in SCORE_COLUMNS.items() if col in df.columns}
missing   = set(SCORE_COLUMNS.keys()) - set(available.keys())
if missing:
    print(f"Warning: Missing scoring columns: {', '.join(missing)}. They will be ignored.")

impact_terms = []

if 'alphamissense_score' in available:
    impact_terms.append(df['alphamissense_score'].fillna(0) * available['alphamissense_score'])

if 'gnomad_max_af' in available:
    impact_terms.append((1 / (1 + df['gnomad_max_af'].fillna(1))) * available.get('gnomad_max_af', 0))

for col in ['spliceai_max_score', 'phylop_score', 'phastcons_score',
            'conservation_score', 'tss_window_score', 'regulatory_score', 'model_score']:
    if col in available:
        impact_terms.append(df[col].fillna(0) * available[col])

if not impact_terms:
    print("Warning: No scoring columns → using random order.")
    df['impact_score'] = pd.Series(np.random.rand(len(df)))
else:
    df['impact_score'] = sum(impact_terms)

# === 3. Select Top N Variants ===
selected = df.sort_values('impact_score', ascending=False).head(TOP_N)
selected_list = selected.to_dict(orient='records')
print(f"Selected {len(selected_list)} variants for testing.")

# === 4. LLM Prompt Template (SAFE VERSION) ===
PROTOCOL_PROMPT_TEMPLATE = """
You are a molecular biology expert designing functional assays for rare disease variants.

Variant: {variant_id}
Gene: {gene}
Change: {protein_change}{consequence_part}
Predicted Mechanism: {mechanism}
Domain: {domain}
AlphaMissense: {alpha_missense:.3f}{cadd_part}, AF: {allele_frequency:.2e}

Design a **single, focused, high-information assay** to test functional impact.
Prioritize: minigene (splicing), trafficking (IF), enzymatic activity, or stability (Western/thermal).
Also make sure that the assays are feasible and budget friendly.

Output in Markdown:
- **Assay Type**
- **Cell Line**
- **Construct Design** (brief)
- **Readout**
- **Positive/Negative Controls**
- **Expected WT vs Mutant Difference**

Keep under 300 words. Be specific and feasible for a CRO.
"""

# === 5. Generate Protocol for Each Variant (SAFE) ===
protocols = []
for i, var in enumerate(selected_list):
    print(f"Generating protocol {i+1}/{len(selected_list)}: {var['variant_id']}")

    safe_var = {
        'variant_id': var.get('variant_id', 'UNKNOWN'),
        'gene': var.get('gene', 'UNKNOWN'),
        'protein_change': var.get('protein_change', 'UNKNOWN'),
        'mechanism': var.get('mechanism', 'unknown'),
        'domain': var.get('domain', 'unknown'),
        'alpha_missense': var.get('alphamissense_score', 0.0),   # note column name
        'allele_frequency': var.get('gnomad_max_af', 1.0),
    }

    consequence = var.get('consequence')
    consequence_part = f" ({consequence})" if consequence else ""

    cadd = var.get('cadd')
    cadd_part = f", CADD: {cadd:.1f}" if cadd is not None else ""

    safe_var.update({'consequence_part': consequence_part, 'cadd_part': cadd_part})

    try:
        prompt = PROTOCOL_PROMPT_TEMPLATE.format(**safe_var)
    except KeyError as e:
        print(f"Prompt formatting error for {var['variant_id']}: {e}")
        continue

    # ---- GROQ CALL ----
    response = client.chat.completions.create(
        model=MODEL,
        messages=[{"role": "user", "content": prompt}],
        temperature=TEMPERATURE,
        max_tokens=600
    )
    # -------------------

    protocol_md = response.choices[0].message.content.strip()

    safe_id = safe_var['variant_id'].replace(":", "_").replace(">", "-").replace("/", "_")
    protocol_file = Path(OUTPUT_DIR) / f"{i+1:02d}_{safe_id}.md"
    protocol_file.write_text(f"# Protocol: {safe_var['variant_id']}\n\n{protocol_md}")

    assay_type = "Unknown"
    for line in protocol_md.split('\n'):
        if line.strip().startswith('**Assay Type**'):
            try:
                assay_type = line.split(":**", 1)[1].strip().split("**")[0].strip()
            except:
                pass
            break

    protocols.append({
        "variant_id": safe_var['variant_id'],
        "protein_change": safe_var['protein_change'],
        "impact_score": round(var.get('impact_score', 0), 4),
        "assay_type": assay_type,
        "protocol_file": protocol_file.name
    })

# === 6. Save Campaign Summary ===
summary = f"""
# Variant Testing Campaign
**Gene**: {selected_list[0].get('gene', 'UNKNOWN')}  
**Date**: {datetime.now().strftime('%Y-%m-%d %H:%M')} GMT  
**Source**: {INPUT_PARQUET}  
**Selected**: {len(selected_list)} / {len(df)} variants  
**Budget**: Top {TOP_N} by impact score  

## Selected Variants
"""
for p in protocols:
    summary += f"- `{p['variant_id']}` | {p['protein_change']} | Score: {p['impact_score']} | **{p['assay_type']}**\n"

summary += f"\n## Protocols Generated\n"
for p in protocols:
    summary += f"- [{p['variant_id']}]({OUTPUT_DIR}/{p['protocol_file']})\n"

Path("campaign_summary.md").write_text(summary)
Path("selected_variants.parquet").write_bytes(selected.to_parquet(index=False))

print(f"\nDone! Check:")
print(f"   - campaign_summary.md")
print(f"   - {OUTPUT_DIR}/ (individual protocols)")
print(f"   - selected_variants.parquet")