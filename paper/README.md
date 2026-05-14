# All Code for the RmvLMM Paper

This folder contains all the R scripts used to generate the results, simulations, and real data analyses presented in our paper. 

These scripts provide a comprehensive **biobank-scale multi-trait GWAS workflow** using the **RmvLMM** (Rotated Multivariate Linear Mixed Model) framework.

## 📂 File Overview

### 1. The Core Engine & Toolbox
* **`MoM_REML_iterative_algorithm.r`**  
  **The Core Algorithm.** This script implements the foundational mathematical model of RmvLMM. It contains the highly optimized iterative algorithm combining the Method of Moments (MoM) and Restricted Maximum Likelihood (REML) to fit the rotated multivariate linear mixed models efficiently.
* **`Testmatrix_pcomb.r`**  
  **The Statistical Toolbox.** An auxiliary script providing essential functions for hypothesis testing, P-value calculations, and P-value. It also includes functions for data generation under H1 and fitting extreme tail distributions.

### 2. Simulation Studies
* **`Power_simulation.r`**  
  **Statistical Power Evaluation.** This script fully simulates the power of the RmvLMM method under various genetic architectures and significance levels.

* **`FPR_simulation.r`**  
  **False Positive Rate Evaluation.** This script evaluates the method's FPR. It involves generating null genotype datasets, calculating test statistics, and fitting/comparing multiple theoretical distributions (Generalized Gamma, Gamma, Chi-square, Log-normal, Weibull).

### 3. Real Data Application
* **`Real_data_analysis-AIO.r`**  
  **Biobank-Scale Multi-Trait GWAS Workflow.** An "All-In-One" (AIO) script that meticulously documents our real-world application on human phenotype data.  
  **💡 Note:** This script is recommended as a reference/template for conducting your own biobank-scale multi-trait GWAS with RmvLMM. It covers the complete pipeline:
  * **Data Processing and Quality Control:** Genotype merging, phenotype normalization (Quantile normalization), covariate adjustment, PLINK-based QC.
  * **GWAS Scanning:** Genome-wide association testing using the RmvLMM method.
  * **Post-GWAS Analysis:** MAGMA-based SNP-to-gene mapping and automated PubMed literature search (via `rentrez`) for significant genes.
  * **Visualization:** High-quality publication-ready plots, including QQ plots, Manhattan plots, Bubble plots for pathway enrichment, and Distribution fitting comparisons.

## 🛠 Dependencies

To run the scripts in this folder, you will need R installed along with the following primary packages:

* **Computation & Parallelization:** `data.table`, `Matrix`, `MASS`, `foreach`, `doParallel`
* **Statistics & Distribution Fitting:** `flexsurv`, `fitdistrplus`, `moments`, `survival`
* **Visualization:** `ggplot2`, `patchwork`, `scales`
* **Bioinformatics & Web Mining:** `rentrez`, `xml2`, `dplyr`, `tidyr`, `stringr`
* **External Tools:** Standard [PLINK](https://www.cog-genomics.org/plink/) and [MAGMA](https://ctg.cncr.nl/software/magma) are required for certain QC and annotation steps in the real data workflow.


*** 
*Note: File paths within the scripts (e.g., `./rtMLMM/data/...`) are configured for our specific server environment. If you wish to reproduce the analyses locally, please adjust the `OUTPUT_DIR` and file paths to match your local directory structure.*
