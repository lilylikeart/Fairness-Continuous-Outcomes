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

## 📦 R Package for Fairness Metrics

To support reproducibility and practical use, we provide an **R package** that implements the proposed fairness metrics.

### Features

- Compute **AreaPDF** and **AreaCDF**
- Kernel Density Estimation (KDE)-based distribution estimation
- Bootstrap confidence intervals for fairness metrics
- Flexible support for different group definitions

### Installation

```r
# Install from GitHub
devtools::install_github("your-username/your-repo-name")
