#' The script provides a simple and rapid example using the `RmvLMM`. 
#' The method for generating simulated data is only for demonstrating how to use `RmvLMM`. For formal research, one must either set more reasonable conditions, or analyze real-world data.

# stopifnot(R.version$major >= 4)
# install.packages("remotes")
# remotes::install_github("amss-stat/RmvLMM")

# ============================================================
# 0. Parameter configuration of demonstration data (For demonstration only)
# ============================================================
n <- 1000    # number of individuals
m <- 2000    # number of SNPs
d <- 10      # number of phenotypes
c <- 2       # number of covariates (eg. intercept + sex)

#' For biobank-scale: 
#' The sample size of a group should be greater than 10,000
#' The number of approximately independent SNPs should be greater than 20,000, which means that the total SNP matrix may be even larger.

# ============================================================
# 1. Simulate SNP genotype matrix G (m x n, SNPs in rows)
# ============================================================
maf <- runif(m, 0.05, 0.5)

G <- matrix(
  rbinom(n * m, 2, rep(maf, each = n)),
  nrow = m,
  ncol = n
)
rownames(G) <- paste0("rs_sim_", 1:m)
colnames(G) <- paste0("Indiv_", 1:n)

# ============================================================
# 2. Prepare sample relatedness matrix K (n x n)
# ============================================================
G_scaled <- scale(t(G))
K <- tcrossprod(G_scaled) / m
K <- K + diag(1e-6, n)

rownames(K) <- paste0("Indiv_", 1:n)
colnames(K) <- paste0("Indiv_", 1:n)

# ============================================================
# 3. Simulate covariate matrix X (n x c)
#    includes intercept (column of 1s) + sex (0/1)
# ============================================================
intercept <- rep(1, n)
sex       <- rbinom(n, 1, 0.5)

X <- cbind(intercept, sex)
rownames(X) <- paste0("Indiv_", 1:n)
colnames(X) <- c("Intercept", "Sex")

# ============================================================
# 4. Simulate phenotype matrix Y (n x d)
# ============================================================
library(MASS)

A     <- matrix(rnorm(d * d), d, d)
Sigma <- A %*% t(A) / d
diag(Sigma) <- diag(Sigma) + 0.5

Y <- mvrnorm(n, mu = rep(0, d), Sigma = Sigma)
rownames(Y) <- paste0("Indiv_", 1:n)
colnames(Y) <- paste0("Trait_", 1:d)

# ============================================================
# 5. Simulate approximately independent SNP list (For demonstration only)
# ============================================================
indep_snps <- data.frame(SNP = rownames(G)[seq(1, m, by = 10)])

write.csv(indep_snps, "independent_snps.csv", row.names = FALSE)

# ============================================================
# 6. Function calls following README
# ============================================================
library(RmvLMM)

# --- Basic Multi-trait GWAS ---
results <- run_RmvLMM(
  Y        = Y,
  X        = X,
  K        = K,
  G        = G,
  out_file = "results.rds"
)

# Results
str(results)
head(results$results)
print(round(results$Vb, 4))
print(round(results$Ve, 4))

df      <- results$results
snp_col <- colnames(df)[1]
p_col   <- colnames(df)[2]  
bonf_thres <- 0.05 / nrow(df)
sig_df <- df[df[[p_col]] < bonf_thres, ]
sig_df <- sig_df[order(sig_df[[p_col]]), ]
print(sig_df)
# The demonstration uses the null hypothesis, so it is not expected that any significant SNPs will be screened out.


# --- Biobank-Scale Aggregation ---
group1 <- 1:500
group2 <- 501:1000

# Group 1: prepare its own K
G_t1      <- t(G[, group1])
G_scaled1 <- scale(G_t1)
K1        <- tcrossprod(G_scaled1) / m
K1        <- K1 + diag(1e-6, nrow(K1))

results_part1 <- run_RmvLMM(
  Y        = Y[group1, ],
  X        = X[group1, ],
  K        = K1,
  G        = G[, group1],
  out_file = "results_part1.rds"
)

# Group 2: prepare its own K
G_t2      <- t(G[, group2])
G_scaled2 <- scale(G_t2)
K2        <- tcrossprod(G_scaled2) / m
K2        <- K2 + diag(1e-6, nrow(K2))

results_part2 <- run_RmvLMM(
  Y        = Y[group2, ],
  X        = X[group2, ],
  K        = K2,
  G        = G[, group2],
  out_file = "results_part2.rds"
)

# --- Combine & Calibrate ---
final_results <- bank_RmvLMM(
  rds_files      = c("results_part1.rds", "results_part2.rds"),
  indep_snp_file = "independent_snps.csv",
  out_file       = "final_results.rds"
)

# View(final_results)
