# RmvLMM: A Rotated Multivariate Linear Mixed Model for Dual Large-Scale GWAS

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R version](https://img.shields.io/badge/R-%3E%3D%204.4.0-blue.svg)](https://www.r-project.org/)
[![Package version](https://img.shields.io/badge/version-0.0.0.9000-orange.svg)](https://github.com/amss-stat/RmvLMM)

Rotated multivariate Linear Mixed Model **(RmvLMM)** is a powerful and scalable statistical framework for dual large-scale GWAS, applicable to biobank-scale samples and a large number of phenotypes. Existing multi-trait GWAS methods are often computationally prohibitive for such massive datasets due to memory constraints and processing time. RmvLMM addresses these challenges using an orthogonal rotation framework, achieving both computational efficiency and high detection power.

## Key Features

- **Orthogonal Rotation:** Decorrelates multiple traits to enable efficient multi-trait analysis and achieve high detection power.
- **Fast MoM-REML Iterative Algorithm:** A novel algorithm combining the Method-of-Moments and REML for rapid covariance matrix estimation.
- **Divided-and-Combined Strategy:** Parallelizes analysis across split-sample groups, aggregating signals into a robust omnibus test and calibrating results via generalized gamma distribution approximation.
- **Biobank Scalability:** Efficiently handles hundreds of thousands of individuals and millions of SNPs. It allows users to flexibly partition samples and process SNP matrices in batches based on available computational resources, enabling efficient multi-trait GWAS at the Biobank scale.
- **Modular Design:** `RmvLMM` package follows a modular philosophy by focusing on the core statistical inference, allowing for seamless integration with high-performance I/O tools like `PLINK`.

## Installation

You can install the development version of RmvLMM from GitHub:

```r
# install.packages("remotes")
remotes::install_github("amss-stat/RmvLMM")
```

### Dependencies
Tested with:
`R (>= 4.4.0)`
`data.table (1.17.2)`
`MASS (7.3.65)`
`survival (3.8.3)`
`flexsurv (2.3.2)`

## Quick Start

### 1. Input Data Format
Before running the analysis, ensure your data follows these formats:

*Note: `RmvLMM` is designed to be compatible with standard GWAS pipelines. Upstream tasks such as quality control and sample/SNP splitting are ideally performed using tools like `PLINK` before loading matrices into `R`.*

**Phenotype Matrix (`Y`)**: $N \times D$ matrix of quantitative traits.
| ID | Trait_1 | Trait_2 | ... | Trait_D |
|:---:|:---:|:---:|:---:|:---:|
| Indiv_1 | 1.25 | -0.42 | ... | 0.88 |
| Indiv_2 | -0.12 | 2.15 | ... | -1.04 |
| Indiv_3 | 0.55 | 0.33 | ... | 0.12 |

**Covariate Matrix (`X`)**: $N \times C$ matrix, must include a column of 1s as the intercept.
| ID | Intercept | Age | Sex | PC1 | ... |
|:---:|:---:|:---:|:---:|:---:|:---:|
| Indiv_1 | 1 | 45 | 1 | -0.012 | ... |
| Indiv_2 | 1 | 52 | 0 | 0.045 | ... |
| Indiv_3 | 1 | 38 | 1 | -0.008 | ... |

**Genotype Matrix (`G`)**: $M \times N$ matrix (SNPs in rows, Individuals in columns). SNP IDs must be provided as row names. A numeric vector is also acceptable for a single SNP.
| SNP | Indiv_1 | Indiv_2 | Indiv_3 | ... |
|:---:|:---:|:---:|:---:|:---:|
| rs1001 | 0 | 1 | 0 | ... |
| rs1002 | 2 | 1 | 0 | ... |
| rs1003 | 0 | 0 | 1 | ... |

*Note:  For all input matrices (`Y, X, G`), IDs should be provided as row/column names and not as the first row/column of the numeric data.*

**Sample Relatedness Matrix (`K`)**: $N \times N$ kinship or GRM matrix. $K$ must be positive definite. If not, consider adding small perturbations to the eigenvalues ​​or diagonal elements.

**Independent SNPs List (`independent_snps.csv`)**: A single-column file containing IDs of at least 20,000 approximately independent SNPs (selected via LD pruning or physical distance). The file must contain a header named SNP with the corresponding IDs.
| SNP |
|:---:|
| rs125 |
| rs458 |
| rs992 |

### 2. Basic Multi-trait GWAS
For smaller datasets (e.g., $N < 20,000$), you can process the entire sample at once:

```r
library(RmvLMM)

# Y: Quantitative phenotype matrix (N x D numeric matrix)
# X: Covariate matrix (N x C numeric matrix, includes a column of 1s)
# K: Sample relatedness matrix (N x N numeric matrix)
# G: Genotype matrix (M x N numeric matrix). SNP IDs must be provided as row names.
#    A numeric vector is also acceptable for a single SNP.
results <- run_RmvLMM(Y = Y, X = X, K = K, G = G, out_file = "results.rds")
```

### 3. Biobank-Scale Aggregation
For large cohorts (e.g., UK Biobank), we recommend a "Divided-and-Combined" strategy:

1.  **Divide Samples:** Divide the total sample into multiple groups, each containing about 10,000 individuals.
2.  **Parallel Processing:** For each group, run `run_RmvLMM()` independently and save the results as `_partXX.rds` files.
3.  **Identify Calibration SNPs:** Prepare an ID list of at least 20,000 approximately independent SNPs.
4.  **Combine & Calibrate:** Use `bank_RmvLMM()` to aggregate the results into a final omnibus test.

```r
# Combine results from multiple group files
final_results <- bank_RmvLMM(
  rds_files = c("_part1.rds", "_part2.rds", "_part3.rds", ...),    # obtained from run_RmvLMM() for all groups
  indep_snp_file = "independent_snps.csv",    # ID List of ~20,000+ approximately independent SNPs
  out_file = "final_calibrated_results.rds"
)
```
*For more details, you can browse the demo code in the `example` folder, or view the function comments in `main.r` in the `R` folder.*

## Performance and Computational Tips

### 1. Memory and Efficiency
RmvLMM is highly optimized for large-scale tasks. 
- **Reference Task:** A task with **26 traits, 10,000 individuals and 100,000 SNPs** typically requires **~16GB of RAM**.
- **Real-world Benchmark:** In our study, an analysis of **323,839 individuals, 454,296 SNPs, and 26 traits** was completed in just **13 hours** using a server with 60 cores and 180GB of RAM.
- **Flexibility:** Users can flexibly adjust the sample group size and SNP batching based on available computational resources. To our knowledge, RmvLMM is currently the only method capable of performing exact multi-trait LMM-based GWAS at this scale within a reasonable timeframe.

### 2. Optimization: Pre-computing Eigen-decomposition
The eigendecomposition of the relatedness matrix $K$ is computationally expensive. If you need to call `run_RmvLMM()` multiple times for the same set of individuals (e.g., when processing SNPs in different batches), you can pre-compute the decomposition and pass it via the `K_precomp` parameter to avoid redundant calculations:

```r
# Pre-compute eigen-decomposition once
eK <- eigen(K, symmetric = TRUE)
K_precomp <- list(
  Ak_diag = eK$values,
  Qk = eK$vectors)

# Pass K_precomp to run_RmvLMM
results <- run_RmvLMM(Y = Y, X = X, G = G, K_precomp = K_precomp, out_file = "batch1.rds")
```

## Citation

If you use RmvLMM in your research, please cite:
> Guo, H.<sup>†</sup>, Zhang, Q.<sup>†</sup>, Zheng, X., Quan, Y., & Li, Q.<sup>*</sup> (2026). *A rotated multivariate linear mixed model for dual large-scale genome-wide association study* (Submitted).
> 
> <sup>†</sup> The first two authors should be regarded as Joint First Authors.  
> <sup>*</sup> Corresponding author. E-mail: [liqz@amss.ac.cn](mailto:liqz@amss.ac.cn)
