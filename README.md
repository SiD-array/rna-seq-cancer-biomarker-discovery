# Simple Biomarker Prediction for Disease Classification

A machine learning pipeline for identifying biomarkers from RNA-Seq gene expression data to classify different cancer types.

## 🎯 Project Overview

This project analyzes RNA-Seq expression data from **801 samples** across **5 cancer types** to identify potential biomarkers for disease classification:

- **BRCA** - Breast Invasive Carcinoma
- **KIRC** - Kidney Renal Clear Cell Carcinoma  
- **LUAD** - Lung Adenocarcinoma
- **PRAD** - Prostate Adenocarcinoma
- **COAD** - Colon Adenocarcinoma

## 📁 Project Structure

```
RNA_Seq_Biomarker_Project/
├── data/
│   ├── raw/                    # Raw data files
│   │   ├── data.csv            # Gene expression matrix (801 × 20,531)
│   │   └── labels.csv          # Cancer type labels
│   └── processed/              # Processed data
│       └── cleaned_scaled_data.pkl
├── scripts/
│   ├── 1_clean_and_eda.py      # Data loading, cleaning, EDA
│   ├── 2_feature_selection.py  # PCA, LASSO, biomarker selection
│   └── 3_train_model.py        # Model training and evaluation
├── src/
│   ├── __init__.py
│   ├── data_loader.py          # Data loading utilities
│   ├── feature_tools.py        # Feature selection functions
│   └── model_trainer.py        # Model training class
├── reports/
│   └── figures/                # Generated visualizations
├── requirements.txt
└── README.md
```

## 🚀 Getting Started

### Prerequisites

- Python 3.8+
- pip

### Installation

1. Clone the repository:
```bash
git clone https://github.com/SiD-array/Simple-Biomarker-Prediction-for-Disease-Classification.git
cd Simple-Biomarker-Prediction-for-Disease-Classification
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Run the pipeline:
```bash
# Step 1: Data Loading & Cleanup
python scripts/1_clean_and_eda.py

# Step 2: Feature Selection (coming soon)
python scripts/2_feature_selection.py

# Step 3: Model Training (coming soon)
python scripts/3_train_model.py
```

## 📊 Dataset Summary

| Metric | Value |
|--------|-------|
| Total Samples | 801 |
| Total Genes | 20,531 |
| Cancer Types | 5 |
| Missing Values | 0% |

### Class Distribution

| Cancer Type | Samples | Percentage |
|-------------|---------|------------|
| BRCA | 300 | 37.5% |
| KIRC | 146 | 18.2% |
| LUAD | 141 | 17.6% |
| PRAD | 136 | 17.0% |
| COAD | 78 | 9.7% |

## 📈 Pipeline Steps

### Step 1: Data Loading & Initial Cleanup ✅
- Load RNA-Seq expression data
- Validate data format and quality
- Check for missing values
- Generate initial visualizations

### Step 2: Feature Selection (In Progress)
- Dimensionality reduction with PCA
- Feature selection using LASSO
- Identify candidate biomarkers

### Step 3: Model Training (Upcoming)
- Train classification models
- Hyperparameter tuning
- Model evaluation and comparison

## 📄 License

This project is for educational purposes.

## 👤 Author

- GitHub: [@SiD-array](https://github.com/SiD-array)



