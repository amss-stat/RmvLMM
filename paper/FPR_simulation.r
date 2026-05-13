# ==============================================================================
# step1. Generate genotype data
# ==============================================================================
# --- Loading libraries and functions ---
library(MASS)
library(foreach)
library(doParallel)
library(data.table)
source("./rtMLMM/code/MoM_REML_iterative_algorithm.r") 
source("./rtMLMM/code/Testmatrix_pcomb.r") 

set.seed(123)

# ===================== 1. Parameters =====================
n_snps <- 1e5
n_samples <- 2e5

maf_min <- 0.1
maf_max <- 0.5

snp_block_size <- 2000      # SNP block, save memory
sample_block_size <- 1e4    # Samples blocks, facilitates writing files.

out_dir <- "./rtMLMM/data/real_data/human/"

prefix <- "geno_larger_part"

# ===================== 2. Generate MAF =====================
maf_vec <- runif(n_snps, maf_min, maf_max)
saveRDS(maf_vec, file.path(out_dir, "MAF_vector_200000.rds"))

# ===================== 3. Initialize file header =====================
n_sample_blocks <- ceiling(n_samples / sample_block_size)

make_sample_names <- function(start, end) {
  sprintf("Sample_%06d", start:end)
}

for (b in seq_len(n_sample_blocks)) {
  s <- (b - 1) * sample_block_size + 1
  e <- min(b * sample_block_size, n_samples)
  header <- c("SNP", make_sample_names(s, e))
  fwrite(as.list(header),
         file.path(out_dir, sprintf("%s%02d_200000.csv", prefix, b)),
         col.names = FALSE)
}

# ===================== 4. Generate SNP block =====================
n_snp_blocks <- ceiling(n_snps / snp_block_size)

for (blk in seq_len(n_snp_blocks)) {

  snp_start <- (blk - 1) * snp_block_size + 1
  snp_end <- min(blk * snp_block_size, n_snps)
  idx <- snp_start:snp_end

  maf_blk <- maf_vec[idx]

  message(sprintf("SNP block %d / %d", blk, n_snp_blocks))

  # -------------------------------------------------
  # For each SNP: generate using rbinom on 300,000 samples at once.
  # -------------------------------------------------
  geno_block <- vapply(
    maf_blk,
    function(p) rbinom(n_samples, size = 2, prob = p),
    integer(n_samples)
  )

  geno_block <- t(geno_block)
  rownames(geno_block) <- sprintf("SNP_%06d", idx)

  # ===================== 5. Split samples when writing to a file =====================
  for (b in seq_len(n_sample_blocks)) {

    s <- (b - 1) * sample_block_size + 1
    e <- min(b * sample_block_size, n_samples)

    out_dt <- data.table(
      SNP = rownames(geno_block),
      geno_block[, s:e, drop = FALSE]
    )

    fwrite(
      out_dt,
      file.path(out_dir, sprintf("%s%02d_300000.csv", prefix, b)),
      append = TRUE,
      col.names = FALSE
    )
  }

  rm(geno_block)
  gc()
}

message("Genotype generation complete.")

# ================================
# step2. Pre-decomposition of the K matrice
# ================================

out_dir <- "./rtMLMM/data/real_data/human/"

cat("\n>>> Step1: Decompose it step by step, and save each part separately ...\n")

for (part in 1:n_parts) {
  cat(sprintf("\n>>> [Part %02d] Read and decompose the K matrix ...\n", part))

  # 1) Reading the K matrix
  K_file   <- sprintf("%s/K_%d_%d_part%02d.csv", out_dir, n, d, part)
  K_matrix <- as.matrix(fread(K_file, header = FALSE))

  # 2) Eigenvalue decomposition
  eig <- eigen(K_matrix, symmetric = TRUE)

  # 3) Regularization
  lambda_min <- min(eig$values)
  if (lambda_min <= 0) {
    adj <- abs(lambda_min) + 1e-8
    diag(K_matrix) <- diag(K_matrix) + adj
    eig <- eigen(K_matrix, symmetric = TRUE)
  }

  # 4) Cholesky decomposition
  L_K <- chol(K_matrix)

  # 5) Results
  res <- list(
    Qk_T    = t(eig$vectors),
    Ak_diag = eig$values,
    L_K     = L_K,
    n       = nrow(K_matrix),
    trK     = sum(eig$values)
  )

  # 6) Save the RDS file separately for the current part.
  saveRDS(res,
          file = file.path(out_dir,
                           sprintf("K_PRECOMP_larger_part%02d.rds", part)))
  cat(sprintf(" [Part %02d] Saved: K_PRECOMP_larger_part%02d.rds\n", part, part))
}

cat("\n>>>  %d parts decomposition saved \n", n_parts)
  

# =====================================================================================
# step3. Simulation of type I error (FPR)
# =====================================================================================

# ---------------------------------------------------------------------------
# 1. Main
# ---------------------------------------------------------------------------
library(MASS)
library(foreach)
library(doParallel)
library(data.table)
source("./rtMLMM/code/MoM_REML_iterative_algorithm.r") 
source("./rtMLMM/code/Testmatrix_pcomb.r") 

cat("=== Start simulation of 30 parts ===\n")
start_time <- Sys.time()

n_parts<-30


# ----------- Initialization -----------
out_dir <- "./rtMLMM/data/real_data/human/"
all_part_results <- vector("list", n_parts)
run_one_repeat <- function(part) {
  cat(sprintf("\n=== Start part %02d ===\n", part))
  
  # 1) Read files
  K_PRECOMP_file <- file.path(out_dir, sprintf("K_PRECOMP_large_part%02d.rds", part))
  K_PRECOMP <- readRDS(K_PRECOMP_file)

  pheno_file <- file.path(out_dir, sprintf("pheno_quantile_larger_part%02d.csv", part))
  Y <- fread(pheno_file, header = TRUE) 
  Y <- as.matrix(Y)
  X <- Y[, 2:4]
  Y <- apply(Y, 2, as.numeric)
  Y <- Y[, 5:24]  
  X <- as.matrix(X)
  X <- apply(X, 2, as.numeric)
  X <- cbind(1, X)
  
  geno_file <- file.path(out_dir, sprintf("geno_larger_part%02d.csv", part))
  g_matrix <- fread(geno_file, header = TRUE)
  g_matrix <- as.matrix(g_matrix[, -1]) 
  
  # 2) LMM Pre-compute
  precomputed_data <- prepare_lmm_precomputation(
    K = NULL,
    X = X,
    Y = Y,
    K_precomp = K_PRECOMP    
  )

  # 3) Fit mvLMM
  results <- fit_multivariate_lmm_optimized(
    precomputed_data,
    max_iter = 200,
    tol = 1e-4,
    verbose = 0
  )

  # 4) Prepare data for testing
  test_data <- prepare_hypothesis_test_data(
    precomputed_data,
    Vb_est = results$Vb,
    Ve_est = results$Ve
  )

  # 5) Rotate G
  Qk_t <- t(precomputed_data$Qk)
  G_rot <- tcrossprod(Qk_t, g_matrix)

  # 6) Calculate p-values for all SNPs 
  all_pvals <- calculate_p_values_for_matrix(
    test_data = test_data,
    G_rot = G_rot,
    comb_method = "Top2",
    mode = "ALL"
  )

  # 7) Save
  save_file <- file.path(out_dir, sprintf("sim_results_part%02d_300000_10.rds", part))
  saveRDS(all_pvals, save_file)
  cat(sprintf("saved: %s\n", save_file))
  
  return(invisible(NULL))  
}


# ---------------------------------------------------------------------------
# 2. Conduct all parts
# ---------------------------------------------------------------------------
cat("=== Start simulation of 30 part ===\n")
start_time <- Sys.time()

for (part in 1:n_parts) {
  run_one_repeat(part)
}

end_time <- Sys.time()
cat("=== Complete. Total time:", round(difftime(end_time, start_time, units = "secs"), 2), "seconds ===\n")

# ---------------------------------------------------------------------------
# 3. Merge all parts
# ---------------------------------------------------------------------------
library(data.table)
library(moments)
library(stats)
source("./rtMLMM/code/Testmatrix_pcomb.r") 


n_parts <- 30
SIGNIFICANCE_LEVEL <- c(
  5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4,
  5e-5, 1e-5, 5e-6, 1e-6, 5e-7, 1e-7, 5e-8, 1e-8
)

# -------------------------
# 1) Read 1st part and get the number of SNPs
# -------------------------
file1 <- file.path(out_dir, "sim_results_part01_300000_10.rds") 
part_res <- readRDS(file1)

top2_p <- part_res$Top2
n_snp  <- length(top2_p)

snp_names <- paste0("SNP_", seq_len(n_snp))

pval_dt <- data.table(SNP = snp_names)

pval_dt[, Part1 := top2_p]

cat("Number of SNPs:", n_snp, "\n")

# -------------------------
# 2) Read other parts and merge by columns
# -------------------------
for (part in 2:n_parts) {
  file_i <- file.path(out_dir, sprintf("sim_results_part%02d_50000_10.rds", part)) #_300000_20
  if (!file.exists(file_i)) stop("Missing files: ", file_i)
  
  part_res <- readRDS(file_i)
  top2_p   <- part_res$Top2
  
  if (length(top2_p) != n_snp) {
    stop(sprintf("Part %02d The number of SNPs is inconsistent.", part))
  }
  
  pval_dt[[paste0("Part", part)]] <- top2_p
  cat(sprintf("Have read Part %02d\n", part))
}

part_cols <- paste0("Part", 1:n_parts)

# -------------------------
# 3) Omnibus: T = -2 * sum(log(p1, p2, ..., pn))
# -------------------------
calc_fisher <- function(pvec) {
  pvec <- pvec[!is.na(pvec)]
  if (length(pvec) == 0) return(NA_real_)
  -2 * sum(log(pvec))
}

# Part1 ~ Part30

pval_dt[, T := apply(.SD, 1, calc_fisher),
        .SDcols = part_cols]


cat("Omnibus calculation complete.\n")

# ---------------------------------------------------------------------------
# 4. Distribution approximation 
# ---------------------------------------------------------------------------
set.seed(321)

n_sample <- 10000
if (nrow(pval_dt) < n_sample) n_sample <- nrow(pval_dt)

sample_idx <- sample(seq_len(nrow(pval_dt)), n_sample, replace = FALSE)

T_valid <- pval_dt$T[sample_idx]
T_valid <- T_valid[!is.na(T_valid)]

cat("Number of randomly selected SNPs: ", length(T_valid), "\n")

T_clean <-T_valid # [T_valid >= lower & T_valid <= upper]

# Generalized Gamma
est_gg <- fit.generalgamma.mle(T_clean, inits = NULL)
mu_gg    <- as.numeric(est_gg["mu"])
sigma_gg <- as.numeric(est_gg["sigma"])
Q_gg     <- as.numeric(est_gg["Q"])

# Gamma
fit_gamma <- fitdistr(T_clean, "gamma")
shape_gamma <- as.numeric(fit_gamma$estimate["shape"])
rate_gamma  <- as.numeric(fit_gamma$estimate["rate"])
scale_gamma <- 1 / rate_gamma

# Log-normal
fit_lognorm <- fitdistr(T_clean[T_clean > 0], "lognormal")
mu_lognorm  <- as.numeric(fit_lognorm$estimate["meanlog"])
sigma_lognorm <- as.numeric(fit_lognorm$estimate["sdlog"])

# Chi-square
negloglik_chisq <- function(df, x) {
  -sum(dchisq(x, df = df, log = TRUE))
}
m <- mean(T_clean)
opt <- optimize(negloglik_chisq, interval = c(1e-3, 10 * m), x = T_clean)
df_chisq <- as.numeric(opt$minimum)

# Weibull
fit_weibull <- fitdistr(T_clean, "weibull")

shape_weibull <- as.numeric(fit_weibull$estimate["shape"])
scale_weibull <- as.numeric(fit_weibull$estimate["scale"])

# ---------------------------------------------------------------------------
# 6. Calculate p-values
# ---------------------------------------------------------------------------
# Generalized Gamma
calc_gengamma_pvalue <- function(x, mu, sigma, Q){
  if(is.na(x)) return(NA)
  pgengamma(q = x, mu = mu, sigma = sigma, Q = Q, lower.tail = FALSE)
}
pval_dt[, gengamma_pvalue := sapply(T, calc_gengamma_pvalue,
                                    mu = mu_gg, sigma = sigma_gg, Q = Q_gg)]

# Gamma
pval_dt[, gamma_pvalue := sapply(T, function(x){
  if(is.na(x)) return(NA)
  pgamma(q = x, shape = shape_gamma, scale = scale_gamma, lower.tail = FALSE)
})]

# Log-normal
pval_dt[, lognorm_pvalue := sapply(T, function(x){
  if(is.na(x) || x <= 0) return(NA)
  plnorm(q = x, meanlog = mu_lognorm, sdlog = sigma_lognorm, lower.tail = FALSE)
})]

# Chi-square
pval_dt[, chisq_pvalue := sapply(T, function(x){
  if(is.na(x) || x < 0) return(NA)
  pchisq(q = x, df = df_chisq, lower.tail = FALSE)
})]

# Weibull
pval_dt[, weibull_pvalue := sapply(T, function(x){
  if(is.na(x) || x < 0) return(NA)
  pweibull(q = x, shape = shape_weibull, scale = scale_weibull, lower.tail = FALSE)
})]

cat("The p-values for all distributions have been calculated.\n")

# ---------------------------------------------------------------------------
# 7. Calculate FPR
# ---------------------------------------------------------------------------
alpha_dt <- data.table(SignificanceLevel = SIGNIFICANCE_LEVEL)
alpha_dt[, `:=`(
  TypeI_PearsonIII = NA_real_,
  TypeI_GenGamma = NA_real_,
  TypeI_Gamma = NA_real_,
  TypeI_LogNorm = NA_real_,
  TypeI_ChiSq = NA_real_,
  TypeI_Weibull = NA_real_
)]

# loop calculation
for (i in seq_along(SIGNIFICANCE_LEVEL)) {
  alpha <- SIGNIFICANCE_LEVEL[i]
  
  alpha_dt[i, TypeI_PearsonIII := mean(pval_dt$pearsonIII_pvalue <= alpha, na.rm = TRUE)]
  alpha_dt[i, TypeI_GenGamma := mean(pval_dt$gengamma_pvalue <= alpha, na.rm = TRUE)]
  alpha_dt[i, TypeI_Gamma := mean(pval_dt$gamma_pvalue <= alpha, na.rm = TRUE)]
  alpha_dt[i, TypeI_LogNorm := mean(pval_dt$lognorm_pvalue <= alpha, na.rm = TRUE)]
  alpha_dt[i, TypeI_ChiSq := mean(pval_dt$chisq_pvalue <= alpha, na.rm = TRUE)]
  alpha_dt[i, TypeI_Weibull := mean(pval_dt$weibull_pvalue <= alpha, na.rm = TRUE)]
}

# Results
cat("=== Print Results ===\n")
print(alpha_dt)

# Save SNP p-values
pval_file <- file.path(OUT_DIR, sprintf("SNP_Pval_Fisher_%s_all.csv", PARAM_TAG))
fwrite(pval_dt, pval_file)

# Save FPR
alpha_file <- file.path(OUT_DIR, sprintf("TypeIError_rates_%s_all.csv", PARAM_TAG))
fwrite(alpha_dt, alpha_file)

cat(" All results have been saved\n")
cat("-", pval_file, "\n")
cat("-", alpha_file, "\n")
cat("-", img_path, "\n")


