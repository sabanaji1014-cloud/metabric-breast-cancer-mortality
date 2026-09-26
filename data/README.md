# Data

The analysis uses the **METABRIC** breast cancer dataset (clinical attributes,
mRNA z-scores for 489 genes and mutation status for 173 genes of 1,904 patients),
as published on Kaggle:

**Breast Cancer Gene Expression Profiles (METABRIC)**
<https://www.kaggle.com/datasets/raghadalharbi/breast-cancer-gene-expression-profiles-metabric>

The data file is not included in this repository. To run the analysis:

1. Download the dataset from the link above.
2. Put `METABRIC_RNA_Mutation.csv` in this folder:

```
data/
└── METABRIC_RNA_Mutation.csv
```

The original data come from the METABRIC study (Curtis et al., 2012, *Nature*;
Pereira et al., 2016, *Nature Communications*) and are also available on
[cBioPortal](https://www.cbioportal.org/study/summary?id=brca_metabric).
