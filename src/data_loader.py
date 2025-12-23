"""
Data Loader Module
==================
Functions for loading and transforming RNA-Seq data.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Tuple, Dict, Any


def load_raw_data(data_path: str, labels_path: str) -> Tuple[pd.DataFrame, pd.Series]:
    """
    Load raw RNA-Seq expression data and corresponding labels.
    
    Parameters
    ----------
    data_path : str
        Path to the expression data CSV file (samples × genes)
    labels_path : str
        Path to the labels CSV file
        
    Returns
    -------
    Tuple[pd.DataFrame, pd.Series]
        - Expression data DataFrame with samples as rows and genes as columns
        - Labels Series with sample IDs as index and cancer types as values
    """
    print("Loading expression data...")
    # Load expression data - first column is sample ID (index)
    data = pd.read_csv(data_path, index_col=0)
    print(f"  ✓ Loaded expression matrix: {data.shape[0]} samples × {data.shape[1]} genes")
    
    print("Loading labels...")
    # Load labels - first column is sample ID (index)
    labels_df = pd.read_csv(labels_path, index_col=0)
    labels = labels_df['Class']
    print(f"  ✓ Loaded labels for {len(labels)} samples")
    
    # Verify alignment
    if not data.index.equals(labels.index):
        # Reorder labels to match data index
        common_samples = data.index.intersection(labels.index)
        print(f"  ⚠ Aligning samples: {len(common_samples)} common samples found")
        data = data.loc[common_samples]
        labels = labels.loc[common_samples]
    
    return data, labels


def check_missing_values(data: pd.DataFrame) -> Dict[str, Any]:
    """
    Check for missing values in the expression data.
    
    Parameters
    ----------
    data : pd.DataFrame
        Expression data DataFrame
        
    Returns
    -------
    Dict[str, Any]
        Dictionary containing missing value statistics:
        - total_missing: Total number of missing values
        - missing_percentage: Percentage of missing values
        - genes_with_missing: Number of genes with at least one missing value
        - samples_with_missing: Number of samples with at least one missing value
        - missing_per_gene: Series with missing count per gene (only genes with missing)
        - missing_per_sample: Series with missing count per sample (only samples with missing)
    """
    print("\nChecking for missing values...")
    
    total_cells = data.size
    total_missing = data.isna().sum().sum()
    missing_percentage = (total_missing / total_cells) * 100
    
    # Missing per gene
    missing_per_gene = data.isna().sum()
    genes_with_missing = (missing_per_gene > 0).sum()
    
    # Missing per sample
    missing_per_sample = data.isna().sum(axis=1)
    samples_with_missing = (missing_per_sample > 0).sum()
    
    results = {
        'total_missing': total_missing,
        'total_cells': total_cells,
        'missing_percentage': missing_percentage,
        'genes_with_missing': genes_with_missing,
        'samples_with_missing': samples_with_missing,
        'missing_per_gene': missing_per_gene[missing_per_gene > 0],
        'missing_per_sample': missing_per_sample[missing_per_sample > 0]
    }
    
    # Print summary
    print(f"  Total cells: {total_cells:,}")
    print(f"  Missing values: {total_missing:,} ({missing_percentage:.4f}%)")
    print(f"  Genes with missing values: {genes_with_missing:,}")
    print(f"  Samples with missing values: {samples_with_missing:,}")
    
    return results


def get_data_summary(data: pd.DataFrame, labels: pd.Series) -> Dict[str, Any]:
    """
    Generate a comprehensive summary of the dataset.
    
    Parameters
    ----------
    data : pd.DataFrame
        Expression data DataFrame
    labels : pd.Series
        Labels Series
        
    Returns
    -------
    Dict[str, Any]
        Dictionary containing dataset summary statistics
    """
    print("\n" + "="*60)
    print("DATA SUMMARY")
    print("="*60)
    
    summary = {
        'n_samples': data.shape[0],
        'n_genes': data.shape[1],
        'class_distribution': labels.value_counts().to_dict(),
        'class_percentages': (labels.value_counts(normalize=True) * 100).to_dict(),
        'data_types': data.dtypes.value_counts().to_dict(),
        'memory_usage_mb': data.memory_usage(deep=True).sum() / (1024**2)
    }
    
    print(f"\n📊 Dataset Shape:")
    print(f"   • Samples: {summary['n_samples']}")
    print(f"   • Genes/Features: {summary['n_genes']}")
    print(f"   • Memory Usage: {summary['memory_usage_mb']:.2f} MB")
    
    print(f"\n📋 Class Distribution:")
    for cls, count in sorted(summary['class_distribution'].items()):
        pct = summary['class_percentages'][cls]
        print(f"   • {cls}: {count} samples ({pct:.1f}%)")
    
    print(f"\n🔢 Data Types:")
    for dtype, count in summary['data_types'].items():
        print(f"   • {dtype}: {count} columns")
    
    # Basic statistics
    print(f"\n📈 Expression Value Statistics:")
    print(f"   • Min value: {data.min().min():.4f}")
    print(f"   • Max value: {data.max().max():.4f}")
    print(f"   • Mean (overall): {data.values.mean():.4f}")
    print(f"   • Std (overall): {data.values.std():.4f}")
    
    summary['min_value'] = data.min().min()
    summary['max_value'] = data.max().max()
    summary['mean_value'] = data.values.mean()
    summary['std_value'] = data.values.std()
    
    print("="*60)
    
    return summary


def check_data_format(data: pd.DataFrame) -> Dict[str, Any]:
    """
    Check the format and quality of the data.
    
    Parameters
    ----------
    data : pd.DataFrame
        Expression data DataFrame
        
    Returns
    -------
    Dict[str, Any]
        Dictionary containing format check results
    """
    print("\n🔍 Checking Data Format...")
    
    results = {
        'all_numeric': True,
        'negative_values': False,
        'infinite_values': False,
        'constant_genes': 0,
        'duplicate_genes': 0,
        'duplicate_samples': 0
    }
    
    # Check if all columns are numeric
    non_numeric = data.select_dtypes(exclude=[np.number]).columns.tolist()
    if non_numeric:
        results['all_numeric'] = False
        results['non_numeric_columns'] = non_numeric
        print(f"   ⚠ Non-numeric columns found: {len(non_numeric)}")
    else:
        print(f"   ✓ All columns are numeric")
    
    # Check for negative values
    if (data < 0).any().any():
        results['negative_values'] = True
        neg_count = (data < 0).sum().sum()
        print(f"   ⚠ Negative values found: {neg_count:,}")
    else:
        print(f"   ✓ No negative values")
    
    # Check for infinite values
    if np.isinf(data.select_dtypes(include=[np.number])).any().any():
        results['infinite_values'] = True
        inf_count = np.isinf(data.select_dtypes(include=[np.number])).sum().sum()
        print(f"   ⚠ Infinite values found: {inf_count:,}")
    else:
        print(f"   ✓ No infinite values")
    
    # Check for constant genes (zero variance)
    gene_variance = data.var()
    constant_genes = (gene_variance == 0).sum()
    results['constant_genes'] = constant_genes
    if constant_genes > 0:
        print(f"   ⚠ Constant genes (zero variance): {constant_genes}")
    else:
        print(f"   ✓ No constant genes")
    
    # Check for duplicate genes (columns)
    duplicate_genes = data.columns.duplicated().sum()
    results['duplicate_genes'] = duplicate_genes
    if duplicate_genes > 0:
        print(f"   ⚠ Duplicate gene names: {duplicate_genes}")
    else:
        print(f"   ✓ No duplicate gene names")
    
    # Check for duplicate samples (rows)
    duplicate_samples = data.index.duplicated().sum()
    results['duplicate_samples'] = duplicate_samples
    if duplicate_samples > 0:
        print(f"   ⚠ Duplicate sample IDs: {duplicate_samples}")
    else:
        print(f"   ✓ No duplicate sample IDs")
    
    return results



