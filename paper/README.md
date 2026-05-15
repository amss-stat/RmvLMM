# All Code for the RmvLMM Paper

This folder contains all the R scripts used to generate the results, simulations, and real data analyses presented in our paper. 

These scripts provide a comprehensive **biobank-scale multi-trait GWAS workflow** using the **RmvLMM** (Rotated Multivariate Linear Mixed Model) framework.

## 📂 File Overview

### 1. The Core Engine & Toolbox
* **`MoM_REML_iterative_algorithm.r`**  
  **The Core Algorithm.** The script implements the core estimation algorithm of RmvLMM. It contains the highly optimized iterative algorithm combining the Method of Moments (MoM) and Restricted Maximum Likelihood (REML) to fit the rotated multivariate linear mixed models efficiently.
* **`Testmatrix_pcomb.r`**  
  **The Statistical Toolbox.** An auxiliary script providing efficient functions for hypothesis testing and p-value calculations. It also includes functions for data generation under H1 and fitting extreme tail distributions .

### 2. Simulation Studies
* **`Power_simulation.r`**  
  **Statistical Power Evaluation.** The script simulates the power of the RmvLMM method under various genetic architectures and significance levels.

* **`FPR_simulation.r`**  
  **False Positive Rate Evaluation.** The script evaluates the method's FPR. It involves generating null genotype datasets, calculating test statistics, and fitting/comparing several theoretical distributions (Generalized Gamma, Gamma, Chi-square, Log-normal, Weibull).

### 3. Real Data Application
* **`Real_data_analysis-AIO.r`**  
  **Biobank-Scale Multi-Trait GWAS Workflow.** An "All-In-One" (AIO) script that meticulously documents our real-world application on human phenotype data.  
  **💡 Note:** The script is recommended as a reference/template for conducting your own biobank-scale multi-trait GWAS with RmvLMM. It covers the complete pipeline:
  * **Data Processing and Quality Control:** Data cleaning, phenotype normalization, covariate adjustment, PLINK-based QC and other I/O operations.
  * **GWAS Scanning:** Genome-wide association testing using the RmvLMM method.
  * **Post-GWAS Analysis:** MAGMA-based SNP-to-gene mapping and automated PubMed literature search for significant genes.
  * **Visualization:** QQ plots, Manhattan plots, Bubble plots for pathway enrichment, and Distribution fitting comparisons.

## 🛠 Dependencies

To run the scripts in this folder, you will need `R-4.4.0` installed along with the following primary packages:

* **Computation & Parallelization:** `data.table 1.17.2`, `Matrix 1.17.3`, `MASS 7.3.65`, `foreach 1.5.2`, `doParallel 1.0.17`
* **Statistics & Distribution Fitting:** `flexsurv 2.3.2`, `fitdistrplus 1.2.2`, `moments 0.14.1`, `survival 3.8.3`
* **Visualization:** `ggplot2 3.5.2`, `patchwork 1.3.0`, `scales 1.4.0`
* **Bioinformatics & Web Mining:** `rentrez 1.2.3`, `xml2 1.3.8`, `dplyr 1.1.4`, `tidyr 1.3.1`, `stringr 1.5.1`
* **External Tools:** `PLINK v1.90b6.24 64-bit (6 Jun 2021)` and `MAGMA v1.10 (linux)` are required for certain QC and annotation steps in the real data workflow.
