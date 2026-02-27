#!/usr/bin/env python3
"""
Advanced Publication-Quality MD Simulation Analysis
With Statistical Comparisons and Distribution Analysis
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path
from scipy import stats

# Publication settings
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'axes.linewidth': 1.2,
    'xtick.major.width': 1.2,
    'ytick.major.width': 1.2,
    'figure.dpi': 300,
    'savefig.dpi': 600,
})

COLORS = {
    'apo': '#2E86AB',
    'test': '#A23B72',
    'standard': '#F18F01'
}

LABELS = {
    'apo': 'Apo',
    'test': 'Berberine',
    'standard': 'Acarbose'
}


def read_xvg(filename):
    """Read GROMACS .xvg file"""
    data_lines = []
    metadata = {'title': '', 'xlabel': '', 'ylabel': ''}
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            if line.startswith('@'):
                if 'title' in line.lower():
                    try: metadata['title'] = line.split('"')[1]
                    except: pass
                elif 'xaxis' in line.lower() and 'label' in line.lower():
                    try: metadata['xlabel'] = line.split('"')[1]
                    except: pass
                elif 'yaxis' in line.lower() and 'label' in line.lower():
                    try: metadata['ylabel'] = line.split('"')[1]
                    except: pass
                continue
            
            try:
                values = [float(x) for x in line.split()]
                data_lines.append(values)
            except ValueError:
                continue
    
    return np.array(data_lines), metadata


def plot_with_stats(ax, files_dict, ylabel, title, convert_factor=1, time_ns=True):
    """Generic plotting function with statistics"""
    all_means = []
    
    for key, filepath in files_dict.items():
        if filepath.exists():
            data, _ = read_xvg(filepath)
            if time_ns:
                x = data[:, 0] / 1000  # to ns
            else:
                x = data[:, 0]
            y = data[:, 1] * convert_factor
            
            ax.plot(x, y, label=LABELS[key], color=COLORS[key], 
                   linewidth=1.5, alpha=0.85)
            
            # Mean of last 50%
            mean_y = np.mean(y[len(y)//2:])
            all_means.append((LABELS[key], mean_y))
            ax.axhline(mean_y, color=COLORS[key], linestyle='--', 
                      alpha=0.3, linewidth=1)
    
    ax.set_xlabel('Time (ns)' if time_ns else 'Residue', fontweight='bold')
    ax.set_ylabel(ylabel, fontweight='bold')
    ax.legend(frameon=True, fontsize=9, loc='best')
    ax.grid(True, alpha=0.2, linestyle='--')
    ax.set_title(title, fontweight='bold', loc='left', fontsize=11)
    
    # Add text box with means
    if all_means:
        text_str = '\n'.join([f'{name}: {val:.2f}' for name, val in all_means])
        ax.text(0.02, 0.98, text_str, transform=ax.transAxes,
               fontsize=8, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))


def create_box_plots(ax, param_name, files_dict, ylabel, convert_factor=1):
    """Create box plots for comparison"""
    data_list = []
    labels_list = []
    colors_list = []
    
    for key, filepath in files_dict.items():
        if filepath.exists():
            data, _ = read_xvg(filepath)
            values = data[:, 1] * convert_factor
            # Use equilibrated part (last 50%)
            eq_values = values[len(values)//2:]
            data_list.append(eq_values)
            labels_list.append(LABELS[key])
            colors_list.append(COLORS[key])
    
    if data_list:
        bp = ax.boxplot(data_list, labels=labels_list, patch_artist=True,
                       widths=0.6, showmeans=True,
                       meanprops=dict(marker='D', markerfacecolor='red', 
                                     markersize=6, markeredgecolor='darkred'))
        
        # Color boxes
        for patch, color in zip(bp['boxes'], colors_list):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        # Perform statistical tests
        if len(data_list) >= 2:
            # T-test between test and standard
            if len(data_list) >= 2:
                t_stat, p_val = stats.ttest_ind(data_list[1], data_list[0])
                ax.text(0.5, 0.98, f'p-value: {p_val:.4f}', 
                       transform=ax.transAxes, ha='center', va='top',
                       fontsize=9, bbox=dict(boxstyle='round', 
                       facecolor='yellow' if p_val < 0.05 else 'white', alpha=0.7))
        
        ax.set_ylabel(ylabel, fontweight='bold')
        ax.grid(True, alpha=0.2, axis='y')
        ax.set_title(f'{param_name} Distribution', fontweight='bold', fontsize=10)


def create_full_analysis():
    """Create comprehensive analysis figure"""
    
    # Define file paths
    files = {
        'rmsd': {
            'apo': Path('rmsd_apo.xvg'),
            'test': Path('rmsd_test.xvg'),
            'standard': Path('rmsd_standard.xvg')
        },
        'rmsf': {
            'apo': Path('rmsf_apo.xvg'),
            'test': Path('rmsf_test.xvg'),
            'standard': Path('rmsf_standard.xvg')
        },
        'rg': {
            'apo': Path('rg_apo.xvg'),
            'test': Path('rg_test.xvg'),
            'standard': Path('rg_standard.xvg')
        },
        'sasa': {
            'apo': Path('sasa_apo.xvg'),
            'test': Path('sasa_test.xvg'),
            'standard': Path('sasa_standard.xvg')
        },
        'hbonds': {
            'test': Path('hbonds_test.xvg'),
            'standard': Path('hbonds_standard.xvg')
        }
    }
    
    # Create main figure
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(4, 3, figure=fig, hspace=0.4, wspace=0.35,
                  left=0.06, right=0.96, top=0.94, bottom=0.05)
    
    # Row 1: RMSD time series and box plot
    ax1 = fig.add_subplot(gs[0, :2])
    ax2 = fig.add_subplot(gs[0, 2])
    
    print("Analyzing RMSD...")
    plot_with_stats(ax1, files['rmsd'], 'RMSD (Å)', '(A) Root Mean Square Deviation', 
                   convert_factor=10)
    create_box_plots(ax2, 'RMSD', files['rmsd'], 'RMSD (Å)', convert_factor=10)
    
    # Row 2: RMSF and Rg
    ax3 = fig.add_subplot(gs[1, :2])
    ax4 = fig.add_subplot(gs[1, 2])
    
    print("Analyzing RMSF...")
    plot_with_stats(ax3, files['rmsf'], 'RMSF (Å)', '(B) Root Mean Square Fluctuation', 
                   convert_factor=10, time_ns=False)
    create_box_plots(ax4, 'Rg', files['rg'], 'Rg (Å)', convert_factor=10)
    
    # Row 3: Rg time series and SASA box plot
    ax5 = fig.add_subplot(gs[2, :2])
    ax6 = fig.add_subplot(gs[2, 2])
    
    print("Analyzing Radius of Gyration...")
    plot_with_stats(ax5, files['rg'], 'Rg (Å)', '(C) Radius of Gyration', 
                   convert_factor=10)
    create_box_plots(ax6, 'SASA', files['sasa'], 'SASA (nm²)')
    
    # Row 4: SASA and H-bonds
    ax7 = fig.add_subplot(gs[3, :2])
    ax8 = fig.add_subplot(gs[3, 2])
    
    print("Analyzing SASA...")
    plot_with_stats(ax7, files['sasa'], 'SASA (nm²)', '(D) Solvent Accessible Surface Area')
    
    print("Analyzing Hydrogen Bonds...")
    plot_with_stats(ax8, files['hbonds'], 'H-bonds', '(E) Hydrogen Bonds', time_ns=True)
    
    # Main title
    fig.suptitle('Comprehensive MD Simulation Analysis: α-Glucosidase with Berberine and Acarbose',
                 fontsize=15, fontweight='bold', y=0.98)
    
    # Save
    output_png = 'MD_Complete_Analysis.png'
    output_svg = 'MD_Complete_Analysis.svg'
    
    plt.savefig(output_png, dpi=600, bbox_inches='tight', facecolor='white')
    plt.savefig(output_svg, format='svg', bbox_inches='tight', facecolor='white')
    
    print(f"\n✓ Saved: {output_png}")
    print(f"✓ Saved: {output_svg}")
    
    plt.close()
    
    return output_png, output_svg


def create_summary_statistics():
    """Create a separate statistics summary figure"""
    
    files_dict = {
        'RMSD (Å)': {'apo': Path('rmsd_apo.xvg'), 'test': Path('rmsd_test.xvg'), 
                      'standard': Path('rmsd_standard.xvg'), 'factor': 10},
        'Rg (Å)': {'apo': Path('rg_apo.xvg'), 'test': Path('rg_test.xvg'), 
                    'standard': Path('rg_standard.xvg'), 'factor': 10},
        'SASA (nm²)': {'apo': Path('sasa_apo.xvg'), 'test': Path('sasa_test.xvg'), 
                        'standard': Path('sasa_standard.xvg'), 'factor': 1},
    }
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle('Statistical Comparison of MD Parameters', fontsize=14, fontweight='bold')
    
    for idx, (param_name, file_info) in enumerate(files_dict.items()):
        factor = file_info.pop('factor')
        create_box_plots(axes[idx], param_name, file_info, param_name, convert_factor=factor)
    
    plt.tight_layout()
    output = 'MD_Statistics_Comparison.png'
    plt.savefig(output, dpi=600, bbox_inches='tight', facecolor='white')
    print(f"✓ Saved: {output}")
    plt.close()
    
    return output


if __name__ == '__main__':
    print("\n" + "="*70)
    print("Advanced MD Simulation Analysis")
    print("="*70 + "\n")
    
    png1, svg1 = create_full_analysis()
    print("\nCreating statistical comparisons...")
    stats_fig = create_summary_statistics()
    
    print("\n" + "="*70)
    print("Analysis Complete!")
    print("="*70)
    print("\nGenerated files:")
    print(f"  1. {png1}")
    print(f"  2. {svg1}")
    print(f"  3. {stats_fig}")
    print("\n✓ All figures are publication-ready (600 DPI)")
