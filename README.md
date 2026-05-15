# RmvLMM: A Rotated Multivariate Linear Mixed Model for Dual Large-Scale GWAS

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R version](https://img.shields.io/badge/R-%3E%3D%204.4.0-blue.svg)](https://www.r-project.org/)

**RmvLMM** is a scalable and powerful statistical framework designed for multi-trait Genome-Wide Association Studies (GWAS) on biobank-scale datasets. 

As detailed in our paper, *"A rotated multivariate linear mixed model for dual large-scale genome-wide association study"*, RmvLMM addresses the computational challenges of analyzing multiple phenotypes across hundreds of thousands of individuals, and achieves high detection power in multi-trait GWAS.

## Key Features

- **Orthogonal Rotation:** Decorrelates multiple traits to enable efficient multivariate analysis.
- **Fast Iterative REML:** A novel algorithm for rapid covariance matrix estimation.
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
For smaller datasets (e.g., $N < 20,000$), use `run_RmvLMM` directly:

```r
library(RmvLMM)

# Y: Phenotype matrix (N x D)
# X: Covariate matrix (N x C)
# K: Sample relatedness matrix (N x N)
# G: Genotype matrix with SNP ID (M x N)
results <- run_RmvLMM(Y = Y, X = X, K = K, G = G, out_file = "group1_results.rds")
```

### 2. Biobank-Scale Aggregation
For massive cohorts, analyze split groups separately and then combine:

```r
# Combine results from different groups and calibrate using independent SNPs
final_results <- bank_RmvLMM(
  rds_files = c("part1.rds", "part2.rds", "part3.rds"),    # obtained from run_RmvLMM() for all groups
  indep_snp_file = "independent_snps.csv",
  out_file = "final_calibrated_results.rds"
)
```

## Repository Structure

- `R/`: Contains the core implementation of the RmvLMM algorithm.
- [`paper/`](./paper/): Contains all original R scripts and simulation codes used in the paper (Power, FPR, and Real Data Analysis workflow).

## Citation

If you use RmvLMM in your research, please cite:
> *A rotated multivariate linear mixed model for dual large-scale genome-wide association study* (Submitted to Nature Genetics).
