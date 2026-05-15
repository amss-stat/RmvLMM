# RmvLMM: A Rotated Multivariate Linear Mixed Model for Dual Large-Scale GWAS

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R version](https://img.shields.io/badge/R-%3E%3D%204.4.0-blue.svg)](https://www.r-project.org/)

Rotated multivariate Linear Mixed Model **(RmvLMM)** is a powerful and scalable statistical framework for dual large-scale GWAS, applicable to biobank-scale samples and a large number of phenotypes. Existing multi-trait GWAS methods are not computationally scalable to such massive datasets, as they are severely constrained by prohibitive demands on computational memory and excessive processing time. RmvLMM addresses the computational challenges and achieves high detection power in multi-trait GWAS.

## Key Features

- **Orthogonal Rotation:** Decorrelates multiple traits to enable efficient multi-traIT analysis.
- **Fast MoM-REML Iterative algorithm:** A novel algorithm for rapid covariance matrix estimation.
- **Divided-and-Combined Strategy:** Supports parallelized analysis across split-sample groups, aggregating signals into a robust omnibus test.
- **Biobank Scalability:** Efficiently handles tens of thousands of individuals and millions of SNPs.

## Installation

You can install the development version of RmvLMM from GitHub:

```r
# install.packages("remotes")
remotes::install_github("amss-stat/RmvLMM")
```

### Dependencies
Tested with:
- R (>= 4.4.0)
- data.table (1.17.2)
- MASS (7.3.65)
- survival (3.8.3)
- flexsurv (2.3.2)

## Quick Start

### 1. Basic Multi-Trait GWAS
For small datasets (e.g., $N < 20,000$), use `run_RmvLMM` directly:

```r
library(RmvLMM)

# Y: Phenotype matrix (N x D)
# X: Covariate matrix (N x C)
# K: Sample relatedness matrix (N x N)
# G: Genotype matrix with SNP ID (M x N)
results <- run_RmvLMM(Y = Y, X = X, K = K, G = G, out_file = "group1_results.rds")
```

### 2. Biobank-Scale Aggregation
For Biobank-scale cohorts, analyze split groups separately and then combine:

```r
# Combine results from different groups and calibrate using approximately independent SNPs
final_results <- bank_RmvLMM(
  rds_files = c("part1.rds", "part2.rds", "part3.rds"),    # obtained from run_RmvLMM() for all groups
  indep_snp_file = "independent_snps.csv",    # an ID list of tens of thousands approximately independent SNPs
  out_file = "final_calibrated_results.rds"
)
```

## Citation

If you use RmvLMM in your research, please cite:
> Guo, H. *et al*. *A rotated multivariate linear mixed model for dual large-scale genome-wide association study* (Submitted to Nature Genetics, 2026).
