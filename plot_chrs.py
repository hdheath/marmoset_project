import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import rcParams

# Set publication-quality parameters using a common font.
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['DejaVu Sans']
rcParams['font.size'] = 12
rcParams['axes.linewidth'] = 1.5

# Define the results directory and ensure it exists.
RESULTS_DIR = "./results"
if not os.path.exists(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

# Define chromosome lengths (in base pairs).
chrom_lengths = {
    "chrY_pri": 15000000,   # 15 Mb for chromosome Y
    "chrX_pri": 150000000   # 150 Mb for chromosome X
}

def plot_chromosome(chrom, ax):
    """
    Plot a chromosome with two slim tracks:
      - Palindrome Arms
      - Composite Sequence Features (colored by source)
    The composite track background is first filled with yellow (indicating 'Ancestral')
    so that any non-annotated (blank) regions appear in yellow.
    All axis spines are removed and the x-axis tick labels (distance) are placed on top.
    """
    chrom_key = f'chr{chrom}_pri'
    chrom_length = chrom_lengths[chrom_key]
    
    # Set slim, closely spaced track centers near the top.
    track_heights = {
        'palindromes': 0.85,   # Palindrome Arms
        'composite': 0.95      # Composite features near the top
    }
    track_height = 0.04  # Very slim track height

    # Set axis limits tightly around the tracks.
    ax.set_xlim(0, chrom_length)
    ax.set_ylim(0, 1)
    
    # Remove all axis spines.
    for spine in ax.spines.values():
        spine.set_visible(False)
    
    # Define a custom color palette.
    color_map = {
        'palindrome': '#e377c2',    # Pink for Palindrome Arms
        'PAR': '#1f77b4',           # Blue
        'X_Transposed': '#ff7f0e',   # Orange
        'Ampliconic': '#2ca02c',     # Green
        'Alpha_Satellite': '#d62728',# Red
        'Satellite': '#9467bd'       # Purple
    }
    
    # Plot Palindrome Arms.
    try:
        palindrome_df = pd.read_csv(os.path.join(RESULTS_DIR, f"chr{chrom}_palindromes.bed"),
                                    sep='\t', names=['chrom', 'start', 'end'])
        for _, row in palindrome_df.iterrows():
            ax.add_patch(patches.Rectangle(
                (row['start'], track_heights['palindromes'] - track_height/2),
                row['end'] - row['start'], track_height,
                facecolor=color_map['palindrome'], edgecolor='none', alpha=0.7
            ))
    except FileNotFoundError as e:
        print(f"Warning: {e.filename} not found for {chrom} palindromes")
    
    # --- Composite Sequence Features Track ---
    # Fill the entire composite track with yellow (for 'Ancestral' regions)
    ax.add_patch(patches.Rectangle(
        (0, track_heights['composite'] - track_height/2),
        chrom_length, track_height,
        facecolor='yellow', edgecolor='none', alpha=0.9
    ))
    
    # Then, overlay composite features from file.
    try:
        composite_df = pd.read_csv(os.path.join(RESULTS_DIR, f"composite_chr{chrom}.bed"),
                                   sep='\t', names=['chrom', 'start', 'end', 'source'])
        for _, row in composite_df.iterrows():
            src = row['source']
            ax.add_patch(patches.Rectangle(
                (row['start'], track_heights['composite'] - track_height/2),
                row['end'] - row['start'], track_height,
                facecolor=color_map.get(src, 'grey'), edgecolor='none', alpha=0.9
            ))
    except FileNotFoundError as e:
        print(f"Warning: {e.filename} not found for {chrom} composite features")
    
    # Set y-axis tick labels (unrotated).
    ax.set_yticks([track_heights['palindromes'], track_heights['composite']])
    ax.set_yticklabels(["Palindrome Arms", "Sequence Features"], va='center')
    
    # Configure x-axis ticks based on chromosome.
    if chrom == 'X':
        xticks = range(0, 150000001, 20000000)
    else:
        xticks = range(0, 15000001, 5000000)
    ax.set_xticks(xticks)
    ax.set_xticklabels([f"{x/1e6:.0f} Mb" for x in xticks])
    
    # Move x-axis tick labels (distance) to the top.
    ax.xaxis.tick_top()
    
    # Add a subtle x-axis grid.
    ax.grid(axis='x', color='#eeeeee', linestyle='-', linewidth=0.8)
    
    return color_map

# -------------------- Plot and Save Chromosome Y --------------------
fig_y, ax_y = plt.subplots(figsize=(15, 8), dpi=100)
color_map_y = plot_chromosome('Y', ax_y)
ax_y.set_title("Y Chromosome Genomic Features", pad=20, fontsize=14, fontweight='bold')

# Build legend for chromosome Y.
legend_elements_y = [
    patches.Patch(facecolor=color_map_y['palindrome'], label='Palindrome Arms'),
    patches.Patch(facecolor='yellow', label='Ancestral')
]
composite_sources = ['PAR', 'X_Transposed', 'Ampliconic', 'Alpha_Satellite', 'Satellite']
for src in composite_sources:
    legend_elements_y.append(patches.Patch(facecolor=color_map_y[src], label=src))

ax_y.legend(handles=legend_elements_y, loc='upper left', bbox_to_anchor=(1.02, 1),
            frameon=False, ncol=2, fontsize=12)

plt.subplots_adjust(right=0.8)
plt.tight_layout(pad=3.0)
plt.savefig(os.path.join(RESULTS_DIR, 'chromosome_features_plot_Y.pdf'),
            dpi=600, bbox_inches='tight')
plt.close(fig_y)

# -------------------- Plot and Save Chromosome X --------------------
fig_x, ax_x = plt.subplots(figsize=(15, 8), dpi=100)
color_map_x = plot_chromosome('X', ax_x)
ax_x.set_title("X Chromosome Genomic Features", pad=20, fontsize=14, fontweight='bold')

legend_elements_x = [
    patches.Patch(facecolor=color_map_x['palindrome'], label='Palindrome Arms'),
    patches.Patch(facecolor='yellow', label='Ancestral')
]
for src in composite_sources:
    legend_elements_x.append(patches.Patch(facecolor=color_map_x[src], label=src))

ax_x.legend(handles=legend_elements_x, loc='upper left', bbox_to_anchor=(1.02, 1),
            frameon=False, ncol=2, fontsize=12)

plt.subplots_adjust(right=0.8)
plt.tight_layout(pad=3.0)
plt.savefig(os.path.join(RESULTS_DIR, 'chromosome_features_plot_X.pdf'),
            dpi=600, bbox_inches='tight')
plt.close(fig_x)
