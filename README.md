# Predicting Breast Cancer Mortality from Clinical and Genomic Data (METABRIC)

![R](https://img.shields.io/badge/R-4.3-276DC3?logo=r&logoColor=white)
![caret](https://img.shields.io/badge/framework-caret-orange)
![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)

Course project for **Statistical Machine Learning** (B.Sc., 7th semester).

The goal is to predict whether a breast cancer patient **dies of the disease**
using clinical data, gene expression (489 genes) and mutation status (173 genes)
from the [METABRIC](https://www.kaggle.com/datasets/raghadalharbi/breast-cancer-gene-expression-profiles-metabric)
cohort. Penalised logistic regression and several non-linear / ensemble models
are compared on three feature sets (non-genetic, genetic, combined), and the
best model is explained with SHAP values.

📄 **Full analysis with code, outputs and plots:
[`analysis/metabric_mortality.md`](analysis/metabric_mortality.md)**

---

## Project structure

```
metabric-breast-cancer-mortality/
├── analysis/
│   ├── metabric_mortality.Rmd   # source: all code and explanations
│   ├── metabric_mortality.md    # rendered report (readable on GitHub)
│   └── figures/                 # plots produced by the report
├── data/
│   └── README.md                # where to download the dataset
├── results/
│   └── model_comparison.csv     # CV and test metrics of all 21 models
├── install_packages.R           # installs the required R packages
├── LICENSE
└── README.md
```

## Methods

**1. Data preparation**

- Target: *Died of Disease* vs *Living* (patients who died of other causes are
  excluded) → 1,423 patients, 43.7% deaths.
- Leakage control: `patient_id`, `overall_survival` (identical to the target)
  and `overall_survival_months` (follow-up time) are removed.
- Mutation columns are recoded from protein-change labels to 0/1 indicators;
  mutations seen in fewer than 10 training patients are dropped.
- Columns with heavy missingness (`tumor_stage`, `3-gene_classifier_subtype`)
  are removed; other missing values are imputed (median / mode).
- Stratified 80/20 train–test split. Imputation, correlation filtering
  (|r| > 0.8) and scaling are **fitted on the training set only**.
- Exploratory plots and PCA.

**2. Models** — all tuned with the same 5-fold CV (metric: ROC-AUC) and
evaluated once on the same held-out test set:

| Model | Tuned hyper-parameters |
|---|---|
| Lasso logistic regression | λ |
| Ridge logistic regression | λ |
| K-nearest neighbours | k |
| SVM (RBF kernel) | σ, C |
| Decision tree | maximum depth |
| Random forest | maximum depth |
| AdaBoost (`gbm`, exponential loss) | number of trees, tree depth |

Each model is trained on three feature sets: **non-genetic** (24 clinical
features), **genetic** (608 expression + mutation features) and **combined**
(632 features).

**3. Interpretation** — SHAP values for the best model, estimated with
permutation sampling (implemented in base R in the report).

## Results

Best feature set for each algorithm (sorted by mean CV ROC-AUC):

| Algorithm | Features | CV ROC-AUC (mean ± SD) | Test ROC-AUC | Test accuracy | Test F1 |
|---|---|---|---|---|---|
| **Lasso** | Combined | **0.765 ± 0.019** | **0.767** | 0.694 | 0.592 |
| AdaBoost | Combined | 0.760 ± 0.021 | 0.759 | 0.673 | 0.587 |
| Ridge | Non-genetic | 0.757 ± 0.039 | 0.749 | 0.690 | 0.614 |
| SVM (RBF) | Non-genetic | 0.756 ± 0.041 | 0.745 | 0.680 | 0.599 |
| KNN | Non-genetic | 0.749 ± 0.040 | 0.744 | 0.690 | 0.542 |
| Random forest | Non-genetic | 0.749 ± 0.044 | 0.749 | 0.680 | 0.547 |
| Decision tree | Non-genetic | 0.695 ± 0.030 | 0.670 | 0.655 | 0.570 |

All 21 model / feature-set combinations are in
[`results/model_comparison.csv`](results/model_comparison.csv).

<p align="center">
  <img src="analysis/figures/roc-curves-1.png" width="48%">
  <img src="analysis/figures/shap-summary-1.png" width="48%">
</p>

**Key findings**

- The best models reach a ROC-AUC of about **0.76–0.77**. The top six
  algorithms are within one standard deviation of each other, so the sparse and
  interpretable **Lasso** model (41 non-zero coefficients) is chosen.
- **Clinical variables carry most of the signal.** Genetic features alone give
  ROC-AUC ≈ 0.70 and add only a small gain on top of the clinical features.
- SHAP ranks **positive lymph nodes, age, *STAT5A* expression, NPI and type of
  surgery** as the most important features; a *GATA3* mutation and the luminal A
  subtype are linked to lower risk.
- With a 0.5 threshold, recall for the *Died* class is about 0.5; the threshold
  should be tuned if missing high-risk patients is costly.

## How to run

1. Install R (≥ 4.1) and the required packages:
   ```r
   source("install_packages.R")
   ```
2. Download `METABRIC_RNA_Mutation.csv` from
   [Kaggle](https://www.kaggle.com/datasets/raghadalharbi/breast-cancer-gene-expression-profiles-metabric)
   and put it in `data/` (see [`data/README.md`](data/README.md)).
3. Render the report (from the `analysis/` folder, or open the `.Rmd` in
   RStudio and click **Knit**):
   ```r
   rmarkdown::render("analysis/metabric_mortality.Rmd")
   ```
   The full run takes about 30–40 minutes on a laptop (most of the time is
   spent on the SVM, random forest and SHAP steps). A fixed seed (`123`) is used
   for the split and the CV folds.

## Data

METABRIC (Molecular Taxonomy of Breast Cancer International Consortium), as
published on Kaggle. The data are not redistributed in this repository.

- Curtis, C. et al. (2012). The genomic and transcriptomic architecture of 2,000
  breast tumours reveals novel subgroups. *Nature*, 486, 346–352.
- Pereira, B. et al. (2016). The somatic mutation profiles of 2,433 breast
  cancers refine their genomic and transcriptomic landscapes. *Nature
  Communications*, 7, 11479.

## Author

**Saba Naji**

## License

This project is released under the [MIT License](LICENSE).
