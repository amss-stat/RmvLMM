# RmvLMM: A Rotated Multivariate Linear Mixed Model for Dual Large-Scale GWAS

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R version](https://img.shields.io/badge/R-%3E%3D%204.4.0-blue.svg)](https://www.r-project.org/)

**RmvLMM** is a scalable and powerful statistical framework designed for multi-trait Genome-Wide Association Studies (GWAS) on biobank-scale datasets. 

As detailed in our paper, *"A rotated multivariate linear mixed model for dual large-scale genome-wide association study"*, RmvLMM addresses the computational challenges of analyzing multiple phenotypes across hundreds of thousands of individuals.

## Key Features

- **Orthogonal Rotation:** Decorrelates multiple traits to enable efficient multivariate analysis.
- **Fast Iterative REML:** A novel algorithm for rapid covariance matrix estimation.
- **Divided-and-Combined Strategy:** Supports parallelized analysis across split-sample groups, aggregating signals into a robust omnibus test.
- **Biobank Scalability:** Efficiently handles tens of thousands of individuals and millions of SNPs.

## Installation

You can install the development version of RmvLMM from GitHub:

```r
# install.packages("devtools")
devtools::install_github("amss-stat/RmvLMM")
```

### Dependencies
Tested with:
- R (>= 4.4.0)
- data.table (1.17.2)
- MASS (7.3.65)
- flexsurv (2.3.2)
- survival (3.8.3)

## Quick Start

### 1. Basic Multi-Trait GWAS
For smaller datasets (e.g., $N < 20,000$), use `run_RmvLMM` directly:

```r
library(RmvLMM)

# Y: Phenotype matrix (N x D)
# X: Covariate matrix (N x C)
# K: Kinship matrix (N x N)
# G: Genotype matrix (M x N)
results <- run_RmvLMM(Y = Y, X = X, K = K, G = G, out_file = "group1_results.rds")
```

### 2. Biobank-Scale Aggregation
For massive cohorts, analyze split groups separately and then combine:

```r
# Combine results from different groups and calibrate using independent SNPs
final_results <- bank_RmvLMM(
  rds_files = c("part1.rds", "part2.rds", "part3.rds"),
  indep_snp_file = "independent_snps.csv",
  out_file = "final_calibrated_results.rds"
)
```

## Repository Structure

- `R/`: Contains the core implementation of the RmvLMM algorithm.
- [`paper/`](./paper/): Contains all original R scripts and simulation codes used in the paper for reproducibility (Power, FPR, and Real Data Analysis workflow).

## Citation

If you use RmvLMM in your research, please cite:
> *A rotated multivariate linear mixed model for dual large-scale genome-wide association study* (2024).

---

### 给您的最后提示：
1.  **文件格式**：`bank_RmvLMM` 输出的结果列表中包含了 `$results` (data.frame), `$Vb` 和 `$Ve`。
2.  **数据类型提示**：我在文档中强调了 `quantitative traits`，这能避免用户尝试用它跑二分类性状（Case/Control）。
3.  **GitHub 文件夹**：确保您的模拟代码放在 `paper/` 目录下，并保留那个专门介绍原始代码的 `paper/README.md`，这会让您的项目看起来非常专业且井井有条。

祝您的论文顺利发表，工具被广泛引用！
