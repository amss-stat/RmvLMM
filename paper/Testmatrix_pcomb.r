# ============================================================================================================================================================
# Testmatrix_pcomb.r
# ============================================================================================================================================================
# This is an auxiliary script that facilitates hypothesis testing, p-value calculation, power simulation, and false positive rate simulation using the RmvLMM method.
# ============================================================================================================================================================

generate_data_under_H1 <- function(n, d, c, L_K, Vb, Ve, gamma_true, beta_true, p_maf) {
    intercept <- rep(1, n)
    if (c > 0) {
        covariates <- matrix(rnorm(n * c), nrow = n, ncol = c)
        X <- cbind(intercept, covariates)
    } else {
        X <- matrix(intercept, ncol = 1)
    }
    L_K <- L_K
    L_Vb <- t(chol(Vb))
    Z_B <- matrix(rnorm(n * d), n, d)
    B <- L_K %*% Z_B %*% t(L_Vb)
    L_Ve <- t(chol(Ve))
    Z_E <- matrix(rnorm(n * d), n, d)
    E <- Z_E %*% t(L_Ve)
    gamma_matrix <- matrix(0, nrow = c + 1, ncol = d)
    if (c > 0 && !is.null(gamma_true)) gamma_matrix[2, ] <- gamma_true
    Y_base <- (X %*% gamma_matrix) + B + E
    g <- rbinom(n, size = 2, prob = p_maf)
    genetic_effect <- g %*% beta_true
    Y_with_effect <- Y_base + genetic_effect

    return(list(Y = Y_with_effect, X = X, g = g))
}


prepare_hypothesis_test_data <- function(precomputed_data, Vb_est, Ve_est) {
  Y_rot <- precomputed_data$Y_rot
  X_rot <- precomputed_data$X_rot
  Ak_diag <- precomputed_data$Ak_diag
  n <- nrow(Y_rot)
  d <- ncol(Y_rot)
  c <- ncol(X_rot)

  # 1. Eigen decomposition and rotation
  eVe <- eigen(Ve_est, symmetric = TRUE)
  Ae_inv_sqrt_vec <- 1 / sqrt(pmax(eVe$values, 1e-12))
  Qe <- eVe$vectors

  matrix_for_eigh <- crossprod(Qe, Vb_est) %*% Qe * tcrossprod(Ae_inv_sqrt_vec)
  eXi <- eigen(matrix_for_eigh, symmetric = TRUE)
  U <- eXi$vectors
  Xi_diag <- pmax(eXi$values, 1e-12)

  rot_matrix <- sweep(Qe, 2, Ae_inv_sqrt_vec, "*") %*% U
  Y_transformed <- Y_rot %*% rot_matrix

  # 2. Weight matrix calculation
  W_sqrt_inv_matrix <- 1 / sqrt(outer(Ak_diag, Xi_diag, "*") + 1)
  W_sq_matrix <- W_sqrt_inv_matrix^2

  # 3. Stacked Q matrix and Weighted Residuals
  Q_stacked <- matrix(NA, nrow = n, ncol = c * d)
  Ew_matrix <- matrix(NA, nrow = n, ncol = d)
  RSS_H0_vec <- numeric(d)

  for (j in 1:d) {
    W_j <- W_sqrt_inv_matrix[, j]

    # Weighted QR decomposition
    X_j_w <- sweep(X_rot, 1, W_j, "*")
    qr_X_w <- qr(X_j_w)
    Q_j <- qr.Q(qr_X_w)

    # Store weighted Q and residuals
    col_idx <- ((j - 1) * c + 1):(j * c)
    Q_stacked[, col_idx] <- sweep(Q_j, 1, W_j, "*")

    y_j_w <- Y_transformed[, j] * W_j
    resid_j <- qr.resid(qr_X_w, y_j_w)

    Ew_matrix[, j] <- resid_j * W_j
    RSS_H0_vec[j] <- sum(resid_j^2)
  }

  return(list(
    W_sq_matrix = W_sq_matrix,
    Q_stacked = Q_stacked,
    Ew_matrix = Ew_matrix,
    RSS_H0_vec = RSS_H0_vec,
    n = n, d = d, c = c
  ))
}


combine_p_values <- function(p_values, mode = "single", comb_method = "Top2") {
  if (is.vector(p_values)) p_values <- matrix(p_values, nrow = 1)
  p_values[is.na(p_values)] <- 1.0
  n_snp <- nrow(p_values)
  d <- ncol(p_values)
  p_values <- pmin(pmax(p_values, .Machine$double.xmin), 1.0)

  # Cauchy
  if (mode == "ALL" || (mode == "single" && comb_method == "Cauchy")) {
    stat <- rowMeans(tan((0.5 - p_values) * pi))
    res_cauchy <- 0.5 - atan(stat) / pi
    if (mode == "single") {
      return(res_cauchy)
    }
  }

  # Fisher
  if (mode == "ALL" || (mode == "single" && comb_method == "Fisher")) {
    X2 <- -2 * rowSums(log(p_values))
    res_fisher <- pchisq(X2, df = 2 * d, lower.tail = FALSE)
    if (mode == "single") {
      return(res_fisher)
    }
  }

  # Top2 (ACAT)
  if (mode == "ALL" || (mode == "single" && comb_method == "Top2")) {
    if (d < 2) {
      if (mode == "single") stop("Top2 requires at least 2 traits.")
      res_top2 <- rep(NA, n_snp)
    } else {
      # Identify 1st and 2nd smallest P-values
      p_min1 <- do.call(pmin, as.data.frame(p_values))
      min_idx <- max.col(-p_values, ties.method = "first")

      linear_idx <- (1:n_snp) + (min_idx - 1) * as.numeric(n_snp)

      P_temp <- p_values
      P_temp[linear_idx] <- Inf
      p_min2 <- do.call(pmin, as.data.frame(P_temp))
      rm(P_temp)

      # ACAT Calculation
      lp1 <- log(p_min1)
      lp2 <- log(p_min2)
      lc <- lp1 + lp2
      la <- lc / 2
      a_val <- pmin(exp(la), 1 - 1e-16)

      log_p1 <- pbeta(a_val, shape1 = 2, shape2 = d - 1, log.p = TRUE)

      if (d == 2) {
        log_p2 <- rep(-Inf, n_snp)
      } else {
        log_1_minus_a <- log1p(-a_val)
        series_sum <- numeric(n_snp)
        for (k in 1:(d - 2)) {
          series_sum <- series_sum + exp(k * log_1_minus_a) / k
        }
        bracket_val <- -la - series_sum
        is_valid <- bracket_val > 1e-300

        log_p2 <- rep(-Inf, n_snp)
        const_term <- log(d) + log(d - 1)
        if (any(is_valid)) {
          log_p2[is_valid] <- const_term + lc[is_valid] + log(bracket_val[is_valid])
        }
      }

      max_l <- pmax(log_p1, log_p2)
      res_top2 <- exp(max_l + log(exp(log_p1 - max_l) + exp(log_p2 - max_l)))
      res_top2 <- pmin(res_top2, 1.0)
    }
    if (mode == "single") {
      return(res_top2)
    }
  }

  if (mode == "ALL") {
    return(list(Cauchy = res_cauchy, Fisher = res_fisher, Top2 = res_top2))
  }
}


calculate_p_values_for_matrix <- function(test_data, G_rot,
                                          comb_method = "Cauchy",
                                          mode = "single",
                                          block_size = 10000) {
  W_sq_matrix <- test_data$W_sq_matrix
  Q_stacked <- test_data$Q_stacked
  Ew_matrix <- test_data$Ew_matrix
  RSS_H0_vec <- test_data$RSS_H0_vec

  n <- test_data$n
  d <- test_data$d
  c <- test_data$c
  m <- ncol(G_rot)

  df1 <- 1
  df2 <- n - c - 1
  P_matrix <- matrix(NA, nrow = m, ncol = d)

  # Block-wise processing
  start_indices <- seq(1, m, by = block_size)

  for (start_idx in start_indices) {
    end_idx <- min(start_idx + block_size - 1, m)
    idx_seq <- start_idx:end_idx

    # Extract block (drop=FALSE prevents dimension loss)
    G_block <- G_rot[, idx_seq, drop = FALSE]

    # Matrix operations for the block
    V_stacked_block <- crossprod(Q_stacked, G_block)
    Numer_Mat_block <- crossprod(G_block, Ew_matrix)^2

    # Calculate G^2 and Denominator term (Immediate release)
    G_block_sq <- G_block^2
    SS_G_W_block <- crossprod(G_block_sq, W_sq_matrix)
    rm(G_block_sq)

    # Loop over traits
    for (j in 1:d) {
      row_idx_Q <- ((j - 1) * c + 1):(j * c)
      V_j <- V_stacked_block[row_idx_Q, , drop = FALSE]

      ss_H_G_w <- colSums(V_j^2)
      ss_G_w <- SS_G_W_block[, j]
      delta_RSS_num <- Numer_Mat_block[, j]
      rss_h0 <- RSS_H0_vec[j]

      ss_g_ortho <- ss_G_w - ss_H_G_w
      is_stable <- ss_g_ortho > (1e-10 * ss_G_w)

      p_vals_vec <- rep(1.0, length(idx_seq))

      if (any(is_stable)) {
        curr_delta <- delta_RSS_num[is_stable]
        curr_ortho <- ss_g_ortho[is_stable]
        curr_RSS_H1 <- pmax(rss_h0 - (curr_delta / curr_ortho), 1e-12)

        F_stat <- (curr_delta / curr_ortho * df2) / curr_RSS_H1
        p_vals_vec[is_stable] <- pf(F_stat, df1, df2, lower.tail = FALSE)
      }
      P_matrix[idx_seq, j] <- p_vals_vec
    }
  }

  combine_p_values(p_values = P_matrix, mode = mode, comb_method = comb_method)
}


combine_p_values_with_raw_p_mat <- function(p_values, mode = "single", comb_method = "Top2") {
  if (is.vector(p_values)) p_values <- matrix(p_values, nrow = 1)
  p_values[is.na(p_values)] <- 1.0
  n_snp <- nrow(p_values)
  d <- ncol(p_values)
  p_values <- pmin(pmax(p_values, .Machine$double.xmin), 1.0)

  # Cauchy
  if (mode == "ALL" || (mode == "single" && comb_method == "Cauchy")) {
    stat <- rowMeans(tan((0.5 - p_values) * pi))
    res_cauchy <- 0.5 - atan(stat) / pi
    if (mode == "single") {
      return(res_cauchy)
    }
  }

  # Fisher
  if (mode == "ALL" || (mode == "single" && comb_method == "Fisher")) {
    X2 <- -2 * rowSums(log(p_values))
    res_fisher <- pchisq(X2, df = 2 * d, lower.tail = FALSE)
    if (mode == "single") {
      return(res_fisher)
    }
  }

  # Top2 (ACAT)
  if (mode == "ALL" || (mode == "single" && comb_method == "Top2")) {
    if (d < 2) {
      if (mode == "single") stop("Top2 requires at least 2 traits.")
      res_top2 <- rep(NA, n_snp)
    } else {
      # Identify 1st and 2nd smallest P-values
      p_min1 <- do.call(pmin, as.data.frame(p_values))
      min_idx <- max.col(-p_values, ties.method = "first")

      linear_idx <- (1:n_snp) + (min_idx - 1) * as.numeric(n_snp)

      P_temp <- p_values
      P_temp[linear_idx] <- Inf
      p_min2 <- do.call(pmin, as.data.frame(P_temp))
      rm(P_temp)

      # ACAT Calculation
      lp1 <- log(p_min1)
      lp2 <- log(p_min2)
      lc <- lp1 + lp2
      la <- lc / 2
      a_val <- pmin(exp(la), 1 - 1e-16)

      log_p1 <- pbeta(a_val, shape1 = 2, shape2 = d - 1, log.p = TRUE)

      if (d == 2) {
        log_p2 <- rep(-Inf, n_snp)
      } else {
        log_1_minus_a <- log1p(-a_val)
        series_sum <- numeric(n_snp)
        for (k in 1:(d - 2)) {
          series_sum <- series_sum + exp(k * log_1_minus_a) / k
        }
        bracket_val <- -la - series_sum
        is_valid <- bracket_val > 1e-300

        log_p2 <- rep(-Inf, n_snp)
        const_term <- log(d) + log(d - 1)
        if (any(is_valid)) {
          log_p2[is_valid] <- const_term + lc[is_valid] + log(bracket_val[is_valid])
        }
      }

      max_l <- pmax(log_p1, log_p2)
      res_top2 <- exp(max_l + log(exp(log_p1 - max_l) + exp(log_p2 - max_l)))
      res_top2 <- pmin(res_top2, 1.0)
    }
    if (mode == "single") {
      return(res_top2)
    }
  }

  if (mode == "ALL") {
    return(list(Cauchy = res_cauchy, Fisher = res_fisher, Top2 = res_top2, p_values = p_values))
  }
}


calculate_p_values_for_matrix_detailed <- function(test_data, G_rot,
                                          comb_method = "Cauchy",
                                          mode = "single",
                                          block_size = 10000) {
  W_sq_matrix <- test_data$W_sq_matrix
  Q_stacked <- test_data$Q_stacked
  Ew_matrix <- test_data$Ew_matrix
  RSS_H0_vec <- test_data$RSS_H0_vec

  n <- test_data$n
  d <- test_data$d
  c <- test_data$c
  m <- ncol(G_rot)

  df1 <- 1
  df2 <- n - c - 1
  P_matrix <- matrix(NA, nrow = m, ncol = d)

  # Block-wise processing
  start_indices <- seq(1, m, by = block_size)

  for (start_idx in start_indices) {
    end_idx <- min(start_idx + block_size - 1, m)
    idx_seq <- start_idx:end_idx

    # Extract block (drop=FALSE prevents dimension loss)
    G_block <- G_rot[, idx_seq, drop = FALSE]

    # Matrix operations for the block
    V_stacked_block <- crossprod(Q_stacked, G_block)
    Numer_Mat_block <- crossprod(G_block, Ew_matrix)^2

    # Calculate G^2 and Denominator term (Immediate release)
    G_block_sq <- G_block^2
    SS_G_W_block <- crossprod(G_block_sq, W_sq_matrix)
    rm(G_block_sq)

    # Loop over traits
    for (j in 1:d) {
      row_idx_Q <- ((j - 1) * c + 1):(j * c)
      V_j <- V_stacked_block[row_idx_Q, , drop = FALSE]

      ss_H_G_w <- colSums(V_j^2)
      ss_G_w <- SS_G_W_block[, j]
      delta_RSS_num <- Numer_Mat_block[, j]
      rss_h0 <- RSS_H0_vec[j]

      ss_g_ortho <- ss_G_w - ss_H_G_w
      is_stable <- ss_g_ortho > (1e-10 * ss_G_w)

      p_vals_vec <- rep(1.0, length(idx_seq))

      if (any(is_stable)) {
        curr_delta <- delta_RSS_num[is_stable]
        curr_ortho <- ss_g_ortho[is_stable]
        curr_RSS_H1 <- pmax(rss_h0 - (curr_delta / curr_ortho), 1e-12)

        F_stat <- (curr_delta / curr_ortho * df2) / curr_RSS_H1
        p_vals_vec[is_stable] <- pf(F_stat, df1, df2, lower.tail = FALSE)
      }
      P_matrix[idx_seq, j] <- p_vals_vec
    }
  }

  combine_p_values_with_raw_p_mat(p_values = P_matrix, mode = mode, comb_method = comb_method)
}


library(survival)
library(flexsurv)

#' Maximum likelihood estimation of the generalized gamma distribution
fit.generalgamma.mle <- function(data, inits = NULL) {
  
  # 1. Data cleaning
  clean_data <- data[data > 1e-8 & is.finite(data)]
  if(length(clean_data) < 10) stop("Insufficient sample size to perform model fitting.")
  
  # 2. Initialization
  if (is.null(inits)) {
    log_data <- log(clean_data)
    
    init_mu <- mean(log_data)
    init_sigma <- sd(log_data)
    
    n <- length(log_data)
    m2 <- sum((log_data - init_mu)^2) / n
    m22 <- max(m2, 1e-8)
    m3 <- sum((log_data - init_mu)^3) / n
    init_Q <- m3 / (m22^(1.5))
    
    # 3. Range of initial values
    if (init_mu > 20) init_mu <- 20
    if (init_sigma < 0.001) init_sigma <- 0.001
    if (init_sigma > 100) init_sigma <- 100
    if (init_Q > 10) init_Q <- 10  
    if (init_Q < -10) init_Q <- -10
    if (abs(init_Q) < 0.01) {
      init_Q <- 0.01 * ifelse(init_Q >= 0, 1, -1)
    }
    
    inits <- c(init_mu, init_sigma, init_Q)
  }
  
  # 4. Non-gradient methods for MLE
  fit <- flexsurvreg(Surv(clean_data) ~ 1, 
                     dist = "gengamma", 
                     inits = inits,
                     method = "Nelder-Mead",
                     hessian = FALSE,
                     control = list(maxit = 10000, reltol = 1e-8))
  
  est <- fit$res[, "est"]
  names(est) <- c("mu", "sigma", "Q")
  
  return(est)
}

calc_gengamma_pvalue <- function(stat, mu, sigma, Q) {
  pgengamma(q = stat, mu = mu, sigma = sigma, Q = Q, lower.tail = FALSE)
}
