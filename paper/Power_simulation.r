# ==============================================================================
# step1. Environment and parameter configuration
# ==============================================================================
# --- Loading libraries and functions ---
library(MASS)
library(foreach)
library(doParallel)
library(data.table)
source("./rtMLMM/code/MoM_REML_iterative_algorithm.r") 
source("./rtMLMM/code/Testmatrix_pcomb.r")

  d <- 10
  ROU <- 0.4
  V_b <- matrix(0.5, nrow = d, ncol = d)
  diag(V_b) <- 1
  V_e <- outer(1:d, 1:d, function(i, j) ROU^abs(i - j))
  diag(V_e) <- 1 
  
  V_e <- apply(V_e, 2, as.numeric)
  V_b <- apply(V_b, 2, as.numeric)	  
  set.seed(123)
  GAMMA_TRUE <- round(runif(d, min = 0.1, max = 0.5), 2)
  n_parts <- 30
  library(data.table)
  params <- list(
    N_SIMULATIONS = 200,
    SIGNIFICANCE_LEVEL = c(5e-2, 1e-2, 5e-3,10^-2.5, 1e-3,5e-4, 10^-3.5, 1e-4, 5e-5,10^-4.5, 1e-5, 5e-6,10^-5.5, 1e-6, 5e-7,10^-6.5, 1e-7, 5e-8,10^-7.5, 1e-8),
    MAF = 0.5,
    N_INDIVIDUALS = 10000,
    N_TRAITS = 20,
    N_COVARIATES = 1,
    GAMMA_TRUE = GAMMA_TRUE,
    VB_TRUE = V_b,
    VE_TRUE = V_e,
    BETA_TRUE = matrix(c(-0.001, 0.017,  0.041,  0.045,  0.018, -0.023, -0.048, -0.041, -0.015,  0.004), nrow = 1) 
  )


fisher_combine <- function(p_vec) {
  p_vec <- p_vec[!is.na(p_vec) & p_vec > 0 & p_vec <= 1]
  if (length(p_vec) == 0) return(NA_real_)
  
  stat <- -2 * sum(log(p_vec))
  p_comb <- pchisq(stat, df = 2 * length(p_vec), lower.tail = FALSE)
  return(p_comb)
}


start_time <- Sys.time()
out_dir <- "./rtMLMM/data/real_data/human/"

start_time <- Sys.time()
cat(sprintf("--- Start power simulation (N=%d, Using pre-decomposed K) ---\n", params$N_SIMULATIONS))
formatted_beta <- paste(round(params$BETA_TRUE, 4), collapse = ", ")
cat("beta_true: [", formatted_beta, "]\n")


all_part_pvals <- vector("list", params$N_SIMULATIONS)

for (i in seq_len(params$N_SIMULATIONS)) {
  all_part_pvals[[i]] <- vector("list", n_parts)
}

library(doParallel)
library(foreach)
n_cores <- 10
cl <- makeCluster(n_cores)
registerDoParallel(cl)

cat("Parallel backend started with", n_cores, "cores\n")


# ====================================================================================
# step2. Main loop
# ====================================================================================

for (i in 1:params$N_SIMULATIONS ) {

  set.seed(i)

  cat("Starting simulation:", i, "\n")

  sim_top2_vec <- numeric(n_parts)

  # ============ Inner parallel ============
  part_results <- foreach(
    part = seq_len(n_parts),
    .packages = c("MASS", "data.table"),
    .export = c(
      "generate_data_under_H1",
      "prepare_lmm_precomputation",
      "fit_multivariate_lmm_optimized",
      "prepare_hypothesis_test_data",
      "calculate_p_values_for_matrix"
    )
  ) %dopar% {


    K_PRECOMP_file <- file.path(
      out_dir,
      sprintf("K_PRECOMP_large_part%02d.rds", part) 
    )

    K_PRECOMP <- readRDS(K_PRECOMP_file)

    current_data <- generate_data_under_H1(
      n = params$N_INDIVIDUALS,
      d = params$N_TRAITS,
      c = params$N_COVARIATES,
      L_K = K_PRECOMP$L_K,
      Vb  = params$VB_TRUE,
      Ve  = params$VE_TRUE,
      gamma_true = params$GAMMA_TRUE,
      beta_true  = params$BETA_TRUE,
      p_maf      = params$MAF
    )

    g_matrix <- matrix(current_data$g, nrow = 1)

    precomputed_data <- prepare_lmm_precomputation(
      K = NULL,
      X = current_data$X,
      Y = current_data$Y,
      K_precomp = K_PRECOMP
    )

    results <- fit_multivariate_lmm_optimized(
      precomputed_data,
      max_iter = 80,
      tol = 1e-4,
      verbose = 0
    )

    test_data <- prepare_hypothesis_test_data(
      precomputed_data,
      Vb_est = results$Vb,
      Ve_est = results$Ve
    )

    G_rot <- tcrossprod(
      t(precomputed_data$Qk),
      g_matrix
    )

    pvals_list <- calculate_p_values_for_matrix(
      test_data   = test_data,
      G_rot       = G_rot,
      comb_method = "Top2",
      mode        = "ALL"
    )

    list(
      part = part,
      Top2 = pvals_list$Top2,
      pvals_list = pvals_list
    )
  }
  # ============ Parallel termination ============

  for (res in part_results) {
    part <- res$part
    sim_top2_vec[part] <- res$Top2
    all_part_pvals[[i]][[part]] <- res$pvals_list
  }

  combined_p <- pchisq(
    -2 * sum(log(sim_top2_vec)),
    df = 2 * n_parts,
    lower.tail = FALSE
  )

  cat(sprintf(
    "\n>>> Sim %d combine p-value (Fisher) = %.6e\n",
    i, combined_p
  ))

  all_part_pvals[[i]]$combined_p <- combined_p
}

stopCluster(cl)
cat("\n=============================\n")
cat("Simulation complete\n")
cat("=============================\n")

end_time <- Sys.time()
total_time <- end_time - start_time

cat("Total time: ", total_time, "\n")

# ============================
# Save results
# ============================
save_file <- sprintf(
  "./rtMLMM/code/result1/%d_power_list_%dtraits_MAF%.3f.rds",
  params$N_INDIVIDUALS * n_parts,
  params$N_TRAITS,
  params$MAF
)

saveRDS(
  list(
    pvals = all_part_pvals,
    start_time = start_time,
    end_time = end_time,
    total_time = total_time
  ),
  file = save_file
)

cat("\n All p-values ​​have been saved to \n", save_file, "\n")



# ======================================================================================
# step3. Power
# ======================================================================================
res_file <- sprintf(
  "./rtMLMM/code/result1/%d_power_list_%dtraits_MAF%.3f.rds",
  params$N_INDIVIDUALS * n_parts,
  params$N_TRAITS,
  params$MAF
)
sim_res <- readRDS(res_file)

combined_pvec <- sapply(sim_res$pvals, function(x) x$combined_p)

rejection_counts <- sapply(params$SIGNIFICANCE_LEVEL,
                           function(a) sum(combined_pvec < a, na.rm = TRUE))

final_powers <- rejection_counts / params$N_SIMULATIONS

power_table <- data.frame(
  Alpha = params$SIGNIFICANCE_LEVEL,
  Rejections = rejection_counts,
  Power = final_powers
)

cat("Final Power:\n")
print(power_table)
cat("==============================================================\n")

output_power_path <- sprintf(
  "./rtMLMM/code/result1/power_%d_%d_%.2f_Top2.csv",
  params$N_INDIVIDUALS * n_parts,
  params$N_TRAITS,
  params$MAF
)

final_power_matrix <- matrix(final_powers,
                              nrow = length(params$SIGNIFICANCE_LEVEL),
                              dimnames = list(
                                paste0("alpha=", params$SIGNIFICANCE_LEVEL),
                                "Fisher_Top2"
                              ))

write.csv(final_power_matrix, file = output_power_path)

cat("The simulation results regarding Power have been saved to:\n", output_power_path, "\n")

