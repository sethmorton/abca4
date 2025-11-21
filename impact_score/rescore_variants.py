

import pandas as pd
import os
from pathlib import Path

# Get the directory where this script is located
SCRIPT_DIR = Path(__file__).parent.absolute()

# Get the project root directory (parent of impact_score)
PROJECT_ROOT = SCRIPT_DIR.parent

# Construct paths relative to script directory
variants_path = SCRIPT_DIR / "variants_scored.parquet"
# upload a scored variants parquet file here  

variants_df = pd.read_parquet(variants_path)

# Print all available columns to see what we have
print("All available columns:")
print(variants_df.columns.tolist())
print(f"\nDataFrame shape: {variants_df.shape}")
print("\nFirst few rows:")
print(variants_df.head())

# Find all columns that contain "score" in their name (case-insensitive)
score_columns = [col for col in variants_df.columns if 'score' in col.lower()]

print(f"\nFound {len(score_columns)} columns with 'score' in their name:")
print(score_columns)


# Create a dataframe with only the score columns
score_variants_df = variants_df[score_columns].copy()

# Save to pickle file (using .df extension as requested)
output_path = SCRIPT_DIR / "score_variants.df"
directory_path = output_path.parent  # Extracts directory path from full file path

# Create directory if it doesn't exist
directory_path.mkdir(parents=True, exist_ok=True)

score_variants_df.to_pickle(output_path)
print(f"\nSaved {len(score_columns)} score columns to {output_path}")
print(f"Shape of score_variants_df: {score_variants_df.shape}")

#step 1 use a llm to calculate the cost score and time score of each variants 
# add your key here 

from openai import OpenAI

# ←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←
# Paste your Groq key here (starts with gsk_)
client = OpenAI(
    api_key=key,          # ← change this
    base_url="https://api.groq.com/openai/v1"   # Groq endpoint
)
# add you api key here. 

def score_wetlab_experiment(description: str):
    response = client.chat.completions.create(
        model="llama-3.3-70b-versatile",   # Updated to current Groq model
        temperature=0.0,
        max_tokens=100,
        messages=[
            {"role": "system", "content": """You are an expert lab biologist who has run >50 variant functional assays.
Output ONLY these exact three lines, nothing else:

Time: X/10
Cost: Y/10
Reason: one-sentence explanation

Time scale: 1 = <1 week, 10 = >1 month
Cost scale: 1 = <$2k, 10 = >$500k

the given defined variant is experimented using multiplexed wet lab asays, around 3-5 at least and even more if needed"""},
            {"role": "user", "content": f"Experiment: {description}"}
        ]
    )
    content = response.choices[0].message.content.strip()
    print(content)
    
    # Parse the response to extract time and cost scores
    time_score = None
    cost_score = None
    reason = ""
    
    for line in content.split('\n'):
        if line.startswith('Time:'):
            try:
                time_score = float(line.split(':')[1].split('/')[0].strip())
            except:
                pass
        elif line.startswith('Cost:'):
            try:
                cost_score = float(line.split(':')[1].split('/')[0].strip())
            except:
                pass
        elif line.startswith('Reason:'):
            reason = line.split(':', 1)[1].strip()
    
    return time_score, cost_score, reason



#step 2 re-rank the variants based on impact score 

# First, let's identify the key columns for variant description
# Common columns might be: variant_id, hgvs_c, hgvs_p, consequence, etc.
variant_id_cols = [col for col in variants_df.columns if any(x in col.lower() for x in ['variant', 'hgvs', 'id', 'position', 'pos'])]
print(f"\nPotential variant identifier columns: {variant_id_cols}")

# Create a function to build variant description from available columns
def build_variant_description(row):
    """Build a detailed variant description from available columns"""
    parts = []
    
    # Priority columns for variant description
    priority_cols = ['variant_id', 'gene', 'protein_change', 'hgvs_p', 'hgvs_c', 
                     'vep_consequence', 'vep_impact', 'clinical_significance', 
                     'ref', 'alt', 'variant_type', 'domain']
    
    # Add priority columns first
    for col in priority_cols:
        if col in row.index:
            value = row[col]
            if pd.notna(value) and value != '':
                parts.append(f"{col}: {value}")
    
    # If we have parts, join them
    if parts:
        return " | ".join(parts)
    else:
        # Fallback: basic variant info
        return f"Variant at position {row.get('pos', 'unknown')}: {row.get('ref', '?')} > {row.get('alt', '?')}"

# Test with first row
print("\nExample variant description:")
print(build_variant_description(variants_df.iloc[0]))

# Add columns for time and cost scores
variants_df['time_score'] = None
variants_df['cost_score'] = None
variants_df['llm_reason'] = None
variants_df['variant_description'] = None

# Process all variants (or use .head(N) for testing with first N variants)
print("\n" + "="*80)
print("Scoring variants with LLM...")
print("="*80)

for index, row in variants_df.head(10).iterrows():  # Change to variants_df.iterrows() to process all
    exp = build_variant_description(row)
    variants_df.at[index, 'variant_description'] = exp
    
    print(f"\n--- Variant {index} ({index+1}/10) ---")
    print(f"Description: {exp}")
    
    time_score, cost_score, reason = score_wetlab_experiment(exp)
    variants_df.at[index, 'time_score'] = time_score
    variants_df.at[index, 'cost_score'] = cost_score
    variants_df.at[index, 'llm_reason'] = reason
    
# Calculate impact score (lower is better for prioritization)
#all scores here: ['alphamissense_score', 'missense_combined_score', 'spliceai_max_score', 'phylop_score', 'phastcons_score', 'conservation_score', 'tss_window_score', 'regulatory_score', 'model_score']
variants_df['impact_score'] = variants_df['time_score'] + variants_df['cost_score'] + variants_df['alphamissense_score'] + variants_df['missense_combined_score'] + variants_df['spliceai_max_score'] + variants_df['phylop_score'] + variants_df['phastcons_score'] + variants_df['conservation_score'] + variants_df['tss_window_score'] + variants_df['regulatory_score'] + variants_df['model_score']

print("\n" + "="*80)
print("Summary of scored variants:")
print("="*80)
print(variants_df[['variant_id', 'gene', 'protein_change', 'time_score', 'cost_score', 'impact_score']].head(10))

#step 3 save the variants_df to a parquet file

# Create output directory if needed
output_dir = SCRIPT_DIR / "output"
output_dir.mkdir(parents=True, exist_ok=True)

# Save the complete rescored dataframe
output_rescored = output_dir / "variants_scored_rescore.parquet"
variants_df.to_parquet(output_rescored)
print(f"✓ Saved complete rescored variants to: {output_rescored}")

# Also save just the score columns with new LLM scores
score_columns_extended = score_columns + ['time_score', 'cost_score', 'impact_score']
score_variants_extended_df = variants_df[score_columns_extended]
output_scores = output_dir / "score_variants_extended.parquet"
score_variants_extended_df.to_parquet(output_scores)
print(f"✓ Saved extended score columns to: {output_scores}")

print(f"\nProcessing complete! Scored {len(variants_df[variants_df['time_score'].notna()])} variants.")
print(f"Project root: {PROJECT_ROOT}")
print(f"Script directory: {SCRIPT_DIR}")
print(f"Output directory: {output_dir}")

