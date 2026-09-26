# Installs the R packages used in analysis/metabric_mortality.Rmd
packages <- c(
  "readr", "dplyr", "tidyr", "ggplot2",           # data handling and plots
  "caret", "glmnet", "pROC", "kernlab", "rpart",  # models and evaluation
  "ranger", "gbm", "doParallel",
  "rmarkdown", "knitr"                            # rendering the report
)

missing <- setdiff(packages, rownames(installed.packages()))
if (length(missing) > 0) install.packages(missing)
