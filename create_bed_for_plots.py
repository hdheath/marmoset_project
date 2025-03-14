#!/usr/bin/env python3
"""
Create Composite BED Files with Separate Additional Tracks
"""

import os
import pandas as pd

# ---------------------------
# Configuration
# ---------------------------
BASE_DIR = os.getcwd()
RESULTS_DIR = os.path.join(BASE_DIR, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

# Define processing priority (highest first)
PRIORITY_ORDER = ['PAR', 'X_Transposed', 'Ampliconic', 'Alpha_Satellite', 'Satellite']

# Additional track paths
TRACK_PATHS = {
    'Y': {
        'repeat': os.path.join("PAR_annotation", "combined_mask_fastas", "merged_Y.bed"),
        'palindrome': os.path.join("non_PAR_annotation", "palindrome_output", "chrY_palindromes.bed")
    },
    'X': {
        'repeat': os.path.join("PAR_annotation", "combined_mask_fastas", "merged_X.bed"),
        'palindrome': os.path.join("non_PAR_annotation", "palindrome_output", "chrX_palindromes.bed")
    }
}

# Chromosome configurations
CHROM_CONFIG = {
    'Y': {
        'sources': {
            'PAR': [os.path.join("PAR_annotation", "PAR_region", "Y_PAR_region.bed")],
            'Ampliconic': [os.path.join("non_PAR_annotation", "amplicon_output", "chrY_pri_amplicons_final.bed")],
            'Satellite': [os.path.join("non_PAR_annotation", "satellite_output", "chrY_satellites_no_alpha_final_noPAR.bed")],
            'Alpha_Satellite': [os.path.join("non_PAR_annotation", "satellite_output", "chrY_alpha_satellites_final_filtered.bed")],
            'X_Transposed': [os.path.join("non_PAR_annotation", "x_transposed_output", "chrY_filtered_alignments.bed")]
        },
        'output': os.path.join(RESULTS_DIR, "composite_chrY.bed")
    },
    'X': {
        'sources': {
            'PAR': [os.path.join("PAR_annotation", "PAR_region", "X_PAR_region.bed")],
            'Satellite': [os.path.join("non_PAR_annotation", "satellite_output", "chrX_satellites_no_alpha_final_noPAR.bed")],
            'Alpha_Satellite': [os.path.join("non_PAR_annotation", "satellite_output", "chrX_alpha_satellites_final_filtered.bed")],
            'X_Transposed': [os.path.join("non_PAR_annotation", "x_transposed_output", "chrX_filtered_alignments.bed")]
        },
        'output': os.path.join(RESULTS_DIR, "composite_chrX.bed")
    }
}

# ---------------------------
# Data Loading Functions
# ---------------------------
def load_bed(file_path, source):
    """Load a BED file with source column"""
    full_path = os.path.join(BASE_DIR, file_path)
    if not os.path.exists(full_path):
        print(f"Warning: {full_path} not found.")
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'source'])
    
    df = pd.read_csv(full_path, sep='\t', header=None, comment='#',
                     usecols=[0, 1, 2], names=['chrom', 'start', 'end'])
    df['source'] = source
    return df

def load_palindromes(file_path):
    """Load palindrome regions with both arms"""
    full_path = os.path.join(BASE_DIR, file_path)
    if not os.path.exists(full_path):
        print(f"Warning: {full_path} not found.")
        return pd.DataFrame()
    
    df = pd.read_csv(full_path, sep='\t', comment='#', header=None,
                     usecols=[0, 1, 2, 3, 5, 6],
                     names=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'])
    
    arms = pd.concat([
        df[['chrom1', 'start1', 'end1']].rename(columns={'chrom1': 'chrom', 'start1': 'start', 'end1': 'end'}),
        df[['chrom2', 'start2', 'end2']].rename(columns={'chrom2': 'chrom', 'start2': 'start', 'end2': 'end'})
    ])
    return arms

# ---------------------------
# Processing Functions
# ---------------------------
def merge_intervals(df, merge_threshold=5000):
    """Merge intervals within the same source"""
    if df.empty:
        return df
    
    merged = []
    df = df.sort_values(['source', 'start'])
    
    for (source, chrom), group in df.groupby(['source', 'chrom']):
        current = None
        for _, row in group.iterrows():
            if current is None:
                current = row.copy()
            else:
                if row['start'] <= current['end'] + merge_threshold:
                    current['end'] = max(current['end'], row['end'])
                else:
                    merged.append(current)
                    current = row.copy()
        if current is not None:
            merged.append(current)
    
    return pd.DataFrame(merged)

def resolve_priority_overlaps(df):
    """Resolve overlaps according to priority hierarchy"""
    processed = []
    df = df.sort_values(by=['source', 'start'], 
                        key=lambda x: x.map({v: i for i, v in enumerate(PRIORITY_ORDER)})
                        if x.name == 'source' else x)
    
    for _, current_row in df.iterrows():
        current_start = current_row['start']
        current_end = current_row['end']
        overlaps = []
        
        for existing in processed:
            if (current_start < existing['end'] and 
                current_end > existing['start'] and 
                PRIORITY_ORDER.index(current_row['source']) > PRIORITY_ORDER.index(existing['source'])):
                overlaps.append(existing)
        
        for overlap in overlaps:
            if current_start < overlap['start']:
                processed.append({
                    'chrom': current_row['chrom'],
                    'start': current_start,
                    'end': overlap['start'],
                    'source': current_row['source']
                })
            if current_end > overlap['end']:
                current_start = overlap['end']
            else:
                current_start = current_end
                break
        
        if current_start < current_end:
            processed.append({
                'chrom': current_row['chrom'],
                'start': current_start,
                'end': current_end,
                'source': current_row['source']
            })
    
    return pd.DataFrame(processed)

# ---------------------------
# Main Processing
# ---------------------------
def process_chromosome(chrom):
    """Process data for a single chromosome"""
    config = CHROM_CONFIG[chrom]
    composite = pd.DataFrame()
    
    # Process main composite track
    for source, paths in config['sources'].items():
        for path in paths:
            composite = pd.concat([composite, load_bed(path, source)], ignore_index=True)
    
    merged = merge_intervals(composite)
    resolved = resolve_priority_overlaps(merged)
    
    # Filter small satellites only for the X chromosome.
    if chrom == 'X':
        satellite_mask = resolved['source'].eq('Satellite') & (resolved['end'] - resolved['start'] <= 50000)
        small_sats = resolved[satellite_mask]
        resolved = resolved[~satellite_mask]
    else:
        small_sats = pd.DataFrame(columns=resolved.columns)
    
    # Save composite track
    final = resolved.sort_values(['chrom', 'start', 'end'])
    final.to_csv(config['output'], sep='\t', index=False, header=False)
    
    # Save additional tracks
    save_additional_tracks(chrom, small_sats)
    
    print(f"\nChromosome {chrom} Processing Report:")
    print(f"- Composite regions: {len(final)}")
    print(f"- Small satellites filtered: {len(small_sats)}")
    print(f"- Final sources: {', '.join(final['source'].unique())}")

def save_additional_tracks(chrom, small_sats):
    """Save non-composite tracks"""
    # Save repeat regions
    repeat_path = os.path.join(BASE_DIR, TRACK_PATHS[chrom]['repeat'])
    if os.path.exists(repeat_path):
        repeat_df = pd.read_csv(repeat_path, sep='\t', usecols=[0,1,2], 
                              names=['chrom', 'start', 'end'])
        repeat_output = os.path.join(RESULTS_DIR, f"chr{chrom}_repeats.bed")
        repeat_df.sort_values(['chrom', 'start', 'end']).to_csv(
            repeat_output, sep='\t', index=False, header=False
        )
        print(f"- Saved repeats: {len(repeat_df)} regions")
    
    # Save palindrome regions
    palindrome_df = load_palindromes(TRACK_PATHS[chrom]['palindrome'])
    if not palindrome_df.empty:
        palindrome_output = os.path.join(RESULTS_DIR, f"chr{chrom}_palindromes.bed")
        palindrome_df.sort_values(['chrom', 'start', 'end']).to_csv(
            palindrome_output, sep='\t', index=False, header=False
        )
        print(f"- Saved palindromes: {len(palindrome_df)} arms")
    
    # Save small satellites (if any)
    if not small_sats.empty:
        small_sat_output = os.path.join(RESULTS_DIR, f"chr{chrom}_small_satellites.bed")
        small_sats[['chrom', 'start', 'end']].to_csv(
            small_sat_output, sep='\t', index=False, header=False
        )
        print(f"- filtered out : {len(small_sats)} small satellite regions")

# ---------------------------
# Execution
# ---------------------------
if __name__ == "__main__":
    print("Starting genomic feature processing...")
    for chromosome in ['Y', 'X']:
        process_chromosome(chromosome)
    print("\nProcessing complete. Output files saved to:", RESULTS_DIR)
