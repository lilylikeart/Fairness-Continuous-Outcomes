# Fairness Evaluation for Continuous Outcomes

This repository provides code for our research on **fairness evaluation in continuous prediction tasks**, with a focus on educational data. While most fairness research has focused on classification settings, many real-world applications (e.g., test scores, performance prediction) involve **continuous outcomes**, requiring alternative evaluation approaches.

## 📌 Overview

We propose a framework for assessing group fairness in continuous prediction settings by comparing **entire outcome distributions** rather than relying on threshold-based or discretized metrics.

Our approach introduces distribution-based fairness measures that capture disparities between groups more comprehensively and avoid information loss caused by binarization.

## 🔍 Key Contributions

- Adaptation of statistical distance measures for fairness evaluation:
  - **Area Between Probability Density Functions (AreaPDF)**  
    (based on Total Variation Distance)
  - **Area Between Cumulative Distribution Functions (AreaCDF)**  
    (based on Wasserstein Distance)

- A principled approach to:
  - Measure fairness in continuous outcomes
  - Avoid arbitrary thresholding
  - Capture distributional disparities across groups

- Empirical validation in educational prediction tasks

## 📊 Why Continuous Fairness Matters

Most existing fairness metrics assume discrete outcomes (e.g., classification). However, in domains like education:

- Outcomes are inherently continuous (e.g., scores, probabilities)
- Binarization can obscure meaningful disparities
- Distributional differences provide richer insight into inequity

Our framework addresses these challenges directly.

## 📦 fairRegression: R Package for Fairness Metrics

**fairRegression** is an R package for evaluating fairness in regression models using distribution-based fairness metrics. It provides methods to quantify differences between prediction distributions across demographic groups, together with bootstrap confidence intervals for statistical inference.

## Features

- Compute the **AreaCDF** fairness metric based on empirical cumulative distribution functions.
- Compute the **AreaPDF** fairness metric based on kernel density estimation.
- Compute **Conditional AreaPDF** by stratifying on a conditioning variable.
- Estimate confidence intervals using:
  - Percentile bootstrap
  - BCa (Bias-Corrected and Accelerated) bootstrap
  - Studentized (bootstrap-t) intervals
- Compare the fairness of two regression models using bootstrap inference.

---

## Installation

### Install from GitHub

```r
install.packages("remotes")

remotes::install_github("lilylikeart/Fairness-Continuous-Outcomes")
```

### Load the package

```r
library(fairRegression)
```

---

## Quick Start

### Simulated Data

```r
set.seed(1)

df <- data.frame(
  score = c(rnorm(100, 0), rnorm(100, 0.5)),
  group = rep(c(0, 1), each = 100)
)
```

### Compute AreaCDF

```r
result <- compute_AreaCDF(
  data = df,
  target_col = "score",
  group_col = "group"
)

result$AreaCDF
```

### Compute AreaPDF

```r
result <- compute_AreaPDF(
  data = df,
  target_col = "score",
  group_col = "group"
)

result$AreaPDF
```

### Compute Conditional AreaPDF

```r
set.seed(1)

df$age <- runif(nrow(df), 18, 65)

result <- compute_AreaPDF_cond(
  data = df,
  target_col = "score",
  condition_col = "age",
  group_col = "group"
)

result$mean_AreaPDF
```

### Bootstrap Confidence Interval

```r
bootstrap_metric_ci(
  df,
  target = "score",
  group = "group",
  metric = "AreaCDF",
  B = 200
)
```

### Compare Two Models

```r
bootstrap_metric_diff_ci(
  df,
  pred_baseline = "baseline_prediction",
  pred_method = "new_prediction",
  group = "group",
  metric = "AreaCDF",
  B = 200
)
```

---

## Main Functions

| Function | Description |
|----------|-------------|
| `compute_AreaCDF()` | Computes the AreaCDF fairness metric. |
| `compute_AreaPDF()` | Computes the AreaPDF fairness metric. |
| `compute_AreaPDF_cond()` | Computes the conditional AreaPDF fairness metric across bins of a conditioning variable. |
| `bootstrap_metric_ci()` | Computes bootstrap confidence intervals for fairness metrics. |
| `bootstrap_metric_diff_ci()` | Compares two regression models using bootstrap inference. |
| `check_bootstrap_skew()` | Computes diagnostic statistics for bootstrap distributions. |

---

## Documentation

After installing the package, documentation for each function is available through:

```r
?compute_AreaCDF
?compute_AreaPDF
?compute_AreaPDF_cond
?bootstrap_metric_ci
?bootstrap_metric_diff_ci
```

---

## Citation

If you use **fairRegression** in academic research, please cite the associated publication (citation information will be added once available).

---

## License

This package is released under the MIT License.
