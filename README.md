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

n <- 200

df <- data.frame(
  actual = rnorm(n, mean = 50, sd = 10),
  prediction = c(
    rnorm(100, mean = 50, sd = 9),
    rnorm(100, mean = 52, sd = 11)
  ),
  group = rep(c(0, 1), each = 100)
)

# Prediction error
df$error <- df$prediction - df$actual

# Absolute prediction error
df$abs_error <- abs(df$error)
```

### Compute AreaCDF

AreaCDF compares the distribution of prediction errors (or absolute prediction errors) between demographic groups.

```r
result <- compute_AreaCDF(
  data = df,
  target_col = "abs_error",
  group_col = "group"
)

result$AreaCDF
```

### Compute AreaPDF

AreaPDF measures the difference between the probability density functions of prediction errors (or absolute prediction errors) across groups.

```r
result <- compute_AreaPDF(
  data = df,
  target_col = "abs_error",
  group_col = "group"
)

result$AreaPDF
```

### Compute Conditional AreaPDF

Conditional AreaPDF compares prediction errors between groups after conditioning on the ground-truth outcome.

```r
result <- compute_AreaPDF_cond(
  data = df,
  target_col = "abs_error",
  condition_col = "actual",
  group_col = "group"
)

result$mean_AreaPDF
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
