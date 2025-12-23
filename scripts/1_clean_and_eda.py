"""
Script 1: Data Loading, Cleaning, and Initial EDA
==================================================
This script handles:
1. Loading the raw RNA-Seq expression data
2. Checking for missing values and data format issues
3. Generating initial data summaries and visualizations

Author: RNA-Seq Biomarker Project
"""

import sys
from pathlib import Path

# Add src directory to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

from src.data_loader import (
    load_raw_data, 
    check_missing_values, 
    get_data_summary,
    check_data_format
)


def plot_class_distribution(labels: pd.Series, save_path: Path):
    """Plot and save class distribution bar chart."""
    plt.figure(figsize=(10, 6))
    
    # Color palette for cancer types
    colors = {
        'BRCA': '#E74C3C',  # Red - Breast cancer
        'PRAD': '#3498DB',  # Blue - Prostate cancer
        'LUAD': '#2ECC71',  # Green - Lung cancer
        'KIRC': '#9B59B6',  # Purple - Kidney cancer
        'COAD': '#F39C12'   # Orange - Colon cancer
    }
    
    class_counts = labels.value_counts().sort_index()
    bar_colors = [colors.get(cls, '#95A5A6') for cls in class_counts.index]
    
    ax = class_counts.plot(kind='bar', color=bar_colors, edgecolor='black', linewidth=1.2)
    
    plt.title('Class Distribution in RNA-Seq Dataset', fontsize=14, fontweight='bold')
    plt.xlabel('Cancer Type', fontsize=12)
    plt.ylabel('Number of Samples', fontsize=12)
    plt.xticks(rotation=0)
    
    # Add value labels on bars
    for i, (idx, val) in enumerate(class_counts.items()):
        ax.text(i, val + 5, str(val), ha='center', va='bottom', fontweight='bold')
    
    # Add percentage labels
    total = len(labels)
    for i, (idx, val) in enumerate(class_counts.items()):
        pct = (val / total) * 100
        ax.text(i, val/2, f'{pct:.1f}%', ha='center', va='center', 
                color='white', fontweight='bold', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved class distribution plot to {save_path}")


def plot_expression_distribution(data: pd.DataFrame, save_path: Path, n_genes: int = 100):
    """Plot expression value distribution for a sample of genes."""
    plt.figure(figsize=(12, 5))
    
    # Sample genes for visualization
    sampled_genes = np.random.choice(data.columns, min(n_genes, len(data.columns)), replace=False)
    sample_data = data[sampled_genes].values.flatten()
    
    # Remove NaN values
    sample_data = sample_data[~np.isnan(sample_data)]
    
    plt.subplot(1, 2, 1)
    plt.hist(sample_data, bins=50, color='#3498DB', edgecolor='black', alpha=0.7)
    plt.title('Expression Value Distribution', fontsize=12, fontweight='bold')
    plt.xlabel('Expression Value')
    plt.ylabel('Frequency')
    
    plt.subplot(1, 2, 2)
    # Log-transformed (add small constant to avoid log(0))
    log_data = np.log1p(sample_data[sample_data >= 0])
    plt.hist(log_data, bins=50, color='#2ECC71', edgecolor='black', alpha=0.7)
    plt.title('Log-Transformed Expression Distribution', fontsize=12, fontweight='bold')
    plt.xlabel('log(Expression + 1)')
    plt.ylabel('Frequency')
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved expression distribution plot to {save_path}")


def main():
    """Main function to run data loading and initial cleanup."""
    
    print("\n" + "="*70)
    print("  RNA-Seq BIOMARKER PROJECT - Step 1: Data Loading & Initial Cleanup")
    print("="*70)
    print(f"  Run Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*70 + "\n")
    
    # Define paths
    project_root = Path(__file__).parent.parent
    data_path = project_root / "data" / "raw" / "data.csv"
    labels_path = project_root / "data" / "raw" / "labels.csv"
    processed_path = project_root / "data" / "processed"
    figures_path = project_root / "reports" / "figures"
    
    # Create directories if they don't exist
    processed_path.mkdir(parents=True, exist_ok=True)
    figures_path.mkdir(parents=True, exist_ok=True)
    
    # =========================================================================
    # STEP 1: Load the Data
    # =========================================================================
    print("STEP 1: Loading Data")
    print("-" * 40)
    data, labels = load_raw_data(str(data_path), str(labels_path))
    
    # =========================================================================
    # STEP 2: Check Data Format
    # =========================================================================
    print("\nSTEP 2: Checking Data Format")
    print("-" * 40)
    format_results = check_data_format(data)
    
    # =========================================================================
    # STEP 3: Check for Missing Values
    # =========================================================================
    print("\nSTEP 3: Checking for Missing Values")
    print("-" * 40)
    missing_results = check_missing_values(data)
    
    # =========================================================================
    # STEP 4: Generate Data Summary
    # =========================================================================
    print("\nSTEP 4: Generating Data Summary")
    print("-" * 40)
    summary = get_data_summary(data, labels)
    
    # =========================================================================
    # STEP 5: Generate Visualizations
    # =========================================================================
    print("\nSTEP 5: Generating Visualizations")
    print("-" * 40)
    
    # Plot class distribution
    plot_class_distribution(labels, figures_path / "class_distribution.png")
    
    # Plot expression distribution
    plot_expression_distribution(data, figures_path / "expression_distribution.png")
    
    # =========================================================================
    # STEP 6: Save Processed Data
    # =========================================================================
    print("\nSTEP 6: Saving Processed Data")
    print("-" * 40)
    
    # Save as pickle for faster loading in subsequent scripts
    output_file = processed_path / "cleaned_scaled_data.pkl"
    
    # Create a dictionary with all data
    processed_data = {
        'data': data,
        'labels': labels,
        'summary': summary,
        'missing_results': missing_results,
        'format_results': format_results,
        'processing_date': datetime.now().isoformat()
    }
    
    pd.to_pickle(processed_data, output_file)
    print(f"  ✓ Saved processed data to {output_file}")
    
    # =========================================================================
    # Summary Report
    # =========================================================================
    print("\n" + "="*70)
    print("  SUMMARY REPORT")
    print("="*70)
    print(f"""
  Data Loading Status: {'✓ SUCCESS' if data is not None else '✗ FAILED'}
  
  Dataset Overview:
    • Total Samples: {summary['n_samples']}
    • Total Genes: {summary['n_genes']}
    • Cancer Types: {len(summary['class_distribution'])}
    
  Data Quality:
    • Missing Values: {missing_results['total_missing']} ({missing_results['missing_percentage']:.4f}%)
    • All Numeric: {'✓ Yes' if format_results['all_numeric'] else '✗ No'}
    • Negative Values: {'⚠ Yes' if format_results['negative_values'] else '✓ None'}
    • Infinite Values: {'⚠ Yes' if format_results['infinite_values'] else '✓ None'}
    • Constant Genes: {format_results['constant_genes']}
    
  Output Files:
    • Processed Data: {output_file}
    • Class Distribution Plot: {figures_path / 'class_distribution.png'}
    • Expression Distribution Plot: {figures_path / 'expression_distribution.png'}
    """)
    print("="*70)
    print("  Step 1 Complete! Proceed to Step 2: Feature Selection")
    print("="*70 + "\n")
    
    return data, labels, summary


if __name__ == "__main__":
    data, labels, summary = main()



