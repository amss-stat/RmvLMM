# ======================================================================================
# step1. Data acquisition and processing
# ======================================================================================

# ---------------------------------------------------------------------------
# 1. Merge genotype data from all chromosomes
# ---------------------------------------------------------------------------
# [Shell Command] The following steps were performed in Linux terminal:
# MERGE_DIR=./rtMLMM/data/real_data/genotype
# MERGE_LIST=${MERGE_DIR}/merge_list.txt

# rm -f ${MERGE_LIST}

# for chr in {2..22}; do
#   echo "${MERGE_DIR}/chr${chr}" >> ${MERGE_LIST}
# done
# ./rtMLMM/plink \
#   --bfile ./rtMLMM/data/real_data/genotype/chr1 \
#   --merge-list ./rtMLMM/data/real_data/genotype/merge_list.txt \
#   --make-bed \
#   --out ./rtMLMM/data/real_data/human/genome_all

# wc -l ./rtMLMM/data/real_data/human/genome_all.bim
# wc -l ./rtMLMM/data/real_data/human/genome_all.fam


# ---------------------------------------------------------------------------
# 2. Processing tabular data
# ---------------------------------------------------------------------------
library(data.table)

pheno_file  <- "./rtMLMM/data/real_data/human/Blood_Count and Biochemistry.csv"
fam_file    <- "./rtMLMM/data/real_data/human/genome_all.fam"  
output_dir  <- "./rtMLMM/data/real_data/human/"

covariate_cols <- c("eid", "sex", "Townsend", "Age at recruitment")

target_pheno_cols <- c(
  "Alanine aminotransferase",
  "Albumin",
  "Alkaline phosphatase",
  "Apolipoprotein A",
  "Apolipoprotein B",
  "Aspartate aminotransferase",
  "C-reactive protein",
  "Calcium",
  "Cholesterol",
  "Creatinine",
  "Cystatin C",
  "Direct bilirubin",
  "Gamma glutamyltransferase",
  "Glucose",
  "Glycated haemoglobin (HbA1c)",
  "HDL cholesterol",
  "IGF-1",
  "LDL direct",
  "Lipoprotein A",
  "Oestradiol",
  "Phosphate",
  "Rheumatoid factor",
  "SHBG",
  "Testosterone",
  "Total bilirubin",
  "Total protein",
  "Triglycerides",
  "Urate",
  "Urea",
  "Vitamin D"
)

extract_cols <- c(covariate_cols, target_pheno_cols)


col_names <- fread(pheno_file, nrows = 0, header = TRUE)
colnames_clean <- trimws(colnames(col_names))
exist_cols <- intersect(extract_cols, colnames_clean)
cat("Number of columns in the data", length(exist_cols), "Skip non-existent columns\n")

dt_pheno <- fread(
  pheno_file,
  select = exist_cols,
  na.strings = c("", "NA", "N/A")
)
setnames(dt_pheno, trimws(colnames(dt_pheno)))


exist_pheno_cols <- setdiff(exist_cols, covariate_cols)
missing_stats <- data.table(
  var_name = exist_pheno_cols,
  missing_pct = sapply(dt_pheno[, ..exist_pheno_cols], function(x) round(mean(is.na(x)) * 100, 2))
)
cat("\n=== ָMissing rate of specified phenotypes ===\n")
print(missing_stats, row.names = FALSE)

# Screening for phenotypes with a missing rate of less than 20%
keep_pheno_cols <- missing_stats[missing_pct < 20, var_name]
keep_all_cols <- c(covariate_cols, keep_pheno_cols)
dt_pheno_filtered <- dt_pheno[, ..keep_all_cols]
cat("\nNum of columns. Covariates + phenotypes missing under 20%:", length(keep_all_cols), "\n")
cat("Number of high missing phenotypes deleted:", length(exist_pheno_cols) - length(keep_pheno_cols), "\n")

# Delete samples with missing information.
dt_pheno_complete <- dt_pheno_filtered[complete.cases(dt_pheno_filtered)]
cat("Sample size with no missing phenotypes or covariates: ", nrow(dt_pheno_complete), "\n")

# Find the intersection
dt_fam <- fread(fam_file, select = c(1,2), header = FALSE)
colnames(dt_fam) <- c("FID", "IID")
geno_eid <- dt_fam$IID

common_eid <- intersect(dt_pheno_complete$eid, geno_eid)
dt_pheno_common <- dt_pheno_complete[eid %in% common_eid]
cat("Final sample size:", nrow(dt_pheno_common), "\n")

fwrite(dt_pheno_common, paste0(output_dir, "final_pheno_raw_larger.csv"))
sample_list <- data.table(FID = dt_pheno_common$eid, IID = dt_pheno_common$eid)
fwrite(sample_list, paste0(output_dir, "final_samples_larger.txt"), sep = "\t", col.names = FALSE)
cat(paste0(" Phenotypes file", output_dir, "final_pheno_raw_larger.csv\n"))
cat(paste0(" Sample list", output_dir, "final_samples_larger.txt\n"))
cat("Final sample size", nrow(dt_pheno_common), "\n")
cat("Final column", paste(keep_all_cols, collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# 3. Extracting samples
# ---------------------------------------------------------------------------
# [Shell Command] The following steps were performed in Linux terminal:
# BED_PREFIX="./rtMLMM/data/real_data/human/genome_all"
# SAMPLE_LIST="./rtMLMM/data/real_data/human/final_samples_larger.txt"
# OUTPUT_PREFIX="./rtMLMM/data/real_data/human/geno_filtered_larger"

# # PLINK Screening samples
# ./rtMLMM/plink --bfile ${BED_PREFIX} \
#       --keep ${SAMPLE_LIST} \
#       --make-bed \
#       --out ${OUTPUT_PREFIX} \
#       --allow-no-sex   

# echo "=== PLINK completed ==="
# echo "Filtered BE files: ${OUTPUT_PREFIX}.bed"



# ---------------------------------------------------------------------------
# 4. Obtain samples again and process phenotypic data.
# ---------------------------------------------------------------------------
library(data.table)
library(moments)

output_dir  <- "./rtMLMM/data/real_data/human/"
pheno_file  <- "./rtMLMM/data/real_data/human/final_pheno_raw_larger.csv"
sample_list_file <- paste0(output_dir, "final_samples_larger.txt")


target_phenotypes <- c(
  "Alanine aminotransferase","Albumin","Alkaline phosphatase","Apolipoprotein A",
  "Apolipoprotein B","Aspartate aminotransferase","C-reactive protein","Calcium",
  "Cholesterol","Creatinine","Cystatin C","Gamma glutamyltransferase","Glucose",
  "Glycated haemoglobin (HbA1c)","HDL cholesterol","IGF-1","LDL direct",
  "Phosphate","SHBG","Testosterone","Total bilirubin","Total protein",
  "Triglycerides","Urate","Urea","Vitamin D"
)

covariate_cols <- c("eid", "sex", "Townsend", "Age at recruitment")

sample_eid <- fread(sample_list_file)$V2
dt_pheno <- fread(pheno_file, select = c(covariate_cols, target_phenotypes))
dt_pheno <- dt_pheno[eid %in% sample_eid]

cat("Sample siez drawn:", nrow(dt_pheno), "\n")


num_covariates <- c("Townsend","Age at recruitment")
dt_pheno[, (num_covariates) := lapply(.SD, function(x) scale(as.numeric(x))[,1]),
         .SDcols = num_covariates]

dt_pheno <- dt_pheno[complete.cases(dt_pheno[, ..covariate_cols])]
cat("Sample size with complete covariates:", nrow(dt_pheno), "\n")

standardize <- function(x){
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if(is.na(s) || s==0) return(rep(0,length(x)))
  (x - mean(x, na.rm=TRUE))/s
}

log_only <- function(x){
  x <- as.numeric(x)
  x[x<=0] <- NA
  log(x)
}

quantile_norm <- function(x){
  x <- as.numeric(x)
  qnorm((rank(x, na.last="keep") - 0.5) / sum(!is.na(x)))
}

cat("\n Phenotypic missing rate (%):\n")
missing_rate <- sapply(dt_pheno[, ..target_phenotypes], function(x) mean(is.na(x))*100)
print(round(missing_rate, 2))

# RINT 
dt_quant <- copy(dt_pheno)
dt_quant[, (target_phenotypes) := lapply(.SD, quantile_norm), .SDcols=target_phenotypes]

fwrite(dt_quant,   paste0(output_dir, "pheno_quantile_larger.csv"))

cat("\n === Phenotypic processing completed ===\n")
cat("Sample size:", nrow(dt_pheno), "\n")
cat("Number of phenotypes:", length(target_phenotypes), "\n")



# ---------------------------------------------------------------------------
# 5. Data quality control
# ---------------------------------------------------------------------------
# [Shell Command] The following steps were performed in Linux terminal:
# ./rtMLMM/plink --bfile ./rtMLMM/data/real_data/human/geno_filtered_larger --mind 0.05 --make-bed --out ./rtMLMM/data/real_data/human/geno_filtered_larger1  

# ./rtMLMM/plink --bfile ./rtMLMM/data/real_data/human/geno_filtered_larger1 --geno 0.05 --make-bed --out ./rtMLMM/data/real_data/human/geno1_larger

# ./rtMLMM/plink --bfile ./rtMLMM/data/real_data/human/geno1_larger --fill-missing-a2 --make-bed --out ./rtMLMM/data/real_data/human/geno2_larger

# ./rtMLMM/plink --bfile ./rtMLMM/data/real_data/human/geno2_larger --maf 0.01 --hwe 1e-6 --make-bed --out ./rtMLMM/data/real_data/human/geno3_larger


# ./rtMLMM/plink --bfile ./rtMLMM/data/real_data/human/geno3_larger \
#     --indep-pairwise 100 5 0.999999 \
#     --out ./rtMLMM/data/real_data/human/geno3_larger_LDpruned

# ./rtMLMM/plink --bfile ./rtMLMM/data/real_data/human/geno3_larger \
#     --extract ./rtMLMM/data/real_data/human/geno3_larger_LDpruned.prune.in \
#     --make-bed \
#     --out ./rtMLMM/data/real_data/human/geno4_larger

# awk '{print $2, $1, $4}' OFS="\t" ./rtMLMM/data/real_data/human/geno4_larger.bim > ./rtMLMM/data/real_data/human/snp_infoH_geno4_larger.txt  


# ---------------------------------------------------------------------------
# 6. Grouping is required for large samples
# ---------------------------------------------------------------------------
# [Shell Command] The following steps were performed in Linux terminal:
# PLINK=/work/home/acd8ad0oht/rtMLMM/plink 
# PREFIX=/work/home/acd8ad0oht/rtMLMM/data/real_data/human/geno4_larger  
# OUTDIR=/work/home/acd8ad0oht/rtMLMM/data/real_data/human
# CHUNK=10000  
# mkdir -p "$OUTDIR"
# cd "$OUTDIR" || exit 
# echo "=== Step 1: Obtain FID/IID. Raw order fixed ==="
# if [ -f "${PREFIX}.fam" ]; then

#   awk '{print $1, $2}' "${PREFIX}.fam" > all_samples.txt
#   echo "Success! FID/IID extracted from ${PREFIX}.fam to $OUTDIR/all_samples.txt"
# else
#   echo "Error: fam file not found. Please check the path."
#   exit 1  
# fi

# echo "=== Step 2: Preparing to group ==="
# TOTAL_LINES=$(wc -l < all_samples.txt)
# NUM_FULL_CHUNKS=$((TOTAL_LINES / CHUNK))
# REMAINDER=$((TOTAL_LINES % CHUNK))
# echo "Total sample size: $TOTAL_LINES"
# echo "Whole parts (each with $CHUNK): $NUM_FULL_CHUNKS"
# echo "Remained sample: $REMAINDER"

# echo "=== Step 3: Group in order ==="
# split -l $CHUNK all_samples.txt samplelist_part_

# echo "Finished."

# echo "=== Step 4: Generate .raw files ==="
# i=1 
# for f in samplelist_part_*; do
#     OUTFILE="${OUTDIR}/geno4_larger_part${i}"
#     echo "---- Export $i  ----"
#     $PLINK \
#         --bfile "${PREFIX}" \
#         --keep "$f" \
#         --recode A \
#         --out "$OUTFILE"
#     echo ">>> ${OUTFILE}.raw generated"
#     i=$((i+1))
# done

# echo "=== Complete PLINK export, approximately ${CHUNK} samples per batch. ==="

# ---------------------------------------------------------------------------
# 7. Genotype grouping
# ---------------------------------------------------------------------------
library(data.table)

out_dir    <- "./rtMLMM/data/real_data/human/"
snp_info_file <- file.path(out_dir, "snp_infoH_geno4_larger.txt")   # SNP  CHR  POS
snp_info <- fread(snp_info_file, col.names = c("SNP","CHR","POS"))

for (i in 1:33) {
  
  rawfile <- file.path(out_dir, sprintf("geno4_larger_part%d.raw", i))
  message("=== Processing: ", rawfile, " ===")
  geno_raw <- fread(rawfile)
  sample_ids <- as.character(geno_raw$IID) 
  geno_cols <- colnames(geno_raw)[7:ncol(geno_raw)]  # Fisrt 6 cols: FID/IID...
  geno_raw <- as.matrix(geno_raw[, -(1:6), with = FALSE])
  
  snp_sub <- snp_info
  geno_raw  <- as.data.table(t(geno_raw ))
  setnames(geno_raw , sample_ids)  
  geno_raw <- cbind(snp_sub, geno_raw )
  out_file <- file.path(out_dir,
                        sprintf("geno4_larger_part%02d.csv", i))
  fwrite(geno_raw, file = out_file)
  
  message(">>> CSV generated: ", out_file)
  
  rm(geno_raw, snp_sub)
  gc(verbose = FALSE)
}

message("=== All CSV generated ===")


# ---------------------------------------------------------------------------
# 8. Phenotype grouping
# ---------------------------------------------------------------------------
library(data.table)

out_dir <- paste0("./rtMLMM/data/real_data/human")
n <- 323839 
d <- "larger"
n_part <- 33    # number of groups

norm_type <- "quantile"

phen_file <- file.path(out_dir, sprintf("pheno_%s_larger.csv", norm_type)) 
phen_data <- fread(phen_file)

sample_id_col <- "eid"  

for(i in 1:n_part){
  
  geno_file <- file.path(out_dir, sprintf("geno4_larger_part%02d.csv", i))

  geno_colnames <- fread(geno_file, check.names = FALSE, nrows = 0)
  sample_ids <- colnames(geno_colnames)[4:ncol(geno_colnames)]
  
  matched_sample_ids <- intersect(sample_ids, phen_data[[sample_id_col]])
  unmatched_count <- length(sample_ids) - length(matched_sample_ids)
  if (unmatched_count > 0) {
    warning(sprintf("Part %02d with %d samples cannot match in phenotype file", i, unmatched_count))
  }

  phen_filtered <- phen_data[get(sample_id_col) %in% matched_sample_ids, ]
  
  out_phen_file <- file.path(out_dir, sprintf("pheno_%s_larger_part%02d.csv", norm_type, i))
  fwrite(phen_filtered, file = out_phen_file)
 
  cat(sprintf("Part %02d finished: output phenotype file -> %s, matched sample size: %d\n", 
              i, out_phen_file, nrow(phen_filtered)))
}


# ======================================================================================
# step2. Calculation and decomposition of similarity matrices K
# ======================================================================================
library(data.table)
library(Matrix)

n <- 323839        
d <- "larger"
n_parts <- 33      
block_size <- 10000

for (part_idx in 1:33) {

  file_path <- sprintf("./rtMLMM/data/real_data/human/geno4_larger_part%02d.csv"
                        , part_idx)
  cat("=== Processing:", file_path, "===\n")

  total_snps <- nrow(fread(file_path, select = 1))   

  n_samples <- ncol(fread(file_path, nrows = 1)) - 3
  
  cat(total_snps, "\n")
  cat(n_samples, "\n")

  K <- matrix(0, n_samples, n_samples)
  
  # Divide into blocks by row
  start_rows <- seq(1, total_snps, by = block_size)

  for (start_row in start_rows) {
    end_row <- min(start_row + block_size - 1, total_snps)

    snp_block <- fread(file_path, skip = start_row - 1, nrows = end_row - start_row + 1)

    # Extract genotype column (first 3 columns are SNP information).
    X <- as.matrix(snp_block[, 4:ncol(snp_block), with = FALSE])

    X <- t(X) 
    X_scaled <- scale(X, center = TRUE, scale = TRUE)
    K <- K + tcrossprod(X_scaled)

    cat(sprintf(" Finished %02d: SNP row %d-%d\n", part_idx, start_row, end_row))
  }

  K <- K / total_snps

  out_path <- sprintf("./rtMLMM/data/real_data/human/K_larger_part%02d.csv", part_idx)
  write.table(K, out_path, row.names = FALSE, col.names = FALSE, sep = ",")
  cat("=== Saved:", out_path, "===\n\n")
}

out_dir <- "./rtMLMM/data/real_data/human/"
n_parts<-33
cat("\n>>> Step1: Decompose K matrices and save each part ...\n")

for (part in 1:n_parts) {
  cat(sprintf("\n>>> [Part %02d] Read and decompose K ...\n", part))

  K_file   <- sprintf("./rtMLMM/data/real_data/human/K_larger_part%02d.csv", part)
  K_matrix <- as.matrix(fread(K_file, header = FALSE))

  eig <- eigen(K_matrix, symmetric = TRUE)

  lambda_min <- min(eig$values)
  if (lambda_min <= 0) {
    adj <- abs(lambda_min) + 1e-8
    diag(K_matrix) <- diag(K_matrix) + adj    # Regularization
    eig <- eigen(K_matrix, symmetric = TRUE)
  }

  L_K <- chol(K_matrix)

  res <- list(
    Qk    = eig$vectors,
    Ak_diag = eig$values,
    L_K     = L_K,
    n       = nrow(K_matrix),
    trK     = sum(eig$values)
  )

  saveRDS(res,
          file = file.path(out_dir,
                           sprintf("K_PRECOMP_large_part%02d.rds", part)))
  cat(sprintf(" [Part %02d] Saved: K_PRECOMP_large_part%02d.rds\n", part, part))
}

cat("\n>>> %d part EVD results saved. \n", n_parts)


# ======================================================================================
# step3. Main
# ======================================================================================

# ---------------------------------------------------------------------------
# 1. Calculation of p-values ​​before merging groups
# ---------------------------------------------------------------------------
library(MASS)
library(foreach)
library(doParallel)
library(data.table)

source("./rtMLMM/code/Testmatrix_pcomb.r") 
source("./rtMLMM/code/MoM_REML_iterative_algorithm.r") 


norm_type <- "quantile"


for (part in 1:33) {

  cat("\n--- Start ---\n")
  start_time <- Sys.time()

  out_dir <- "./rtMLMM/data/real_data/human"
  
  pheno_file <- file.path(out_dir, sprintf("pheno_%s_larger_part%02d.csv", norm_type, part))
  Y <- fread(pheno_file, header = TRUE) 
  Y <- as.matrix(Y)
  X <- Y[, 2:4]
  Y <- apply(Y, 2, as.numeric)
  Y <- Y[, 5:30]    # remove eid
  X <- as.matrix(X)
  X <- apply(X, 2, as.numeric)
  X <- cbind(1, X)


  params <- list(
      OUTPUT_DIR = "./rtMLMM/code/result/",
      N_INDIVIDUALS = nrow(Y),
      N_TRAITS = ncol(Y),
      N_COVARIATES = ncol(X)-1,
      Y = Y,
      X = X
  )

  batch_size <- 454296    # Number of SNPs processed in batches
  input_file <- sprintf("./rtMLMM/data/real_data/human/geno4_larger_part%02d.csv", part)

  output_file <- sprintf("./rtMLMM/data/real_data/human/result_larger_%s_part%02d.csv", norm_type, part)

  header <- fread(input_file, header = TRUE, nrows = 1)
  total_snp_rows <- as.integer(system(sprintf("wc -l < %s", input_file), intern = TRUE)) - 1
  total_batches <- ceiling(total_snp_rows / batch_size)
  cat(sprintf("Part %02d: Total SNPs = %d. Processing in %d batches, with a maximum of %d SNPs per batch.\n", 
            part, total_snp_rows, total_batches, batch_size))


  # ---------------------------------------
  part_batch_data <- list()

  for (batch_idx in 1:total_batches) {
    skip_rows <- 1 + (batch_idx - 1) * batch_size
    read_rows <- ifelse(batch_idx == total_batches, 
                        total_snp_rows - (batch_idx - 1) * batch_size, 
                        batch_size)
    
    cat(sprintf("Processing Part %02d - Batch %03d: Skipping %d rows, reading %d SNPs...\n", 
            part, batch_idx, skip_rows, read_rows))
    
    g_matrix_batch <- fread(input_file, 
                            header = FALSE, 
                            skip = skip_rows, 
                            nrows = read_rows,
                            col.names = colnames(header))
    
    if (nrow(g_matrix_batch) == 0) {
      warning(sprintf("Part %02d - Batch %03d has no information. Skip.", part, batch_idx))
      next
    }
    
    annot_cols <- g_matrix_batch[, 1:3, with = FALSE]
    g_data <- as.matrix(g_matrix_batch[, 4:ncol(g_matrix_batch), with = FALSE])
    rm(g_matrix_batch)
    gc()
    
    K_PRECOMP_file <- file.path(out_dir, sprintf("K_PRECOMP_large_part%02d.rds", part))
    K_PRECOMP <- readRDS(K_PRECOMP_file)

    precomputed_data <- prepare_lmm_precomputation(
        K = NULL,
        X = params$X,
        Y = params$Y,
        K_precomp = K_PRECOMP
    )
    
    results <- fit_multivariate_lmm_optimized(
        precomputed_data,
        max_iter = 200,
        tol = 1e-4,
        verbose = 0
    )
    
    test_data <- prepare_hypothesis_test_data(
        precomputed_data,
        Vb_est = results$Vb,
        Ve_est = results$Ve
    )
    
  
    Qk_t <- t(precomputed_data$Qk)
    G_rot <- tcrossprod(Qk_t, g_data)
    
    p_value_output2 <- calculate_p_values_for_matrix_detailed(
        test_data = test_data, 
        G_rot = G_rot,
        comb_method = "Top2",
        mode = "ALL"
    )
    
    obj_names <- names(p_value_output2)
    p_name <- obj_names[length(obj_names)]
    method_names <- obj_names[-length(obj_names)]
    
    combined_p <- p_value_output2[method_names]
    
    d <- ncol(p_value_output2[[p_name]])
    p_names <- paste0("p", 1:d)
    
    pval_df <- cbind(
        do.call(cbind, combined_p),
        p_value_output2[[p_name]]
    )
    
    colnames(pval_df) <- c(method_names, p_names)
    final_df <- cbind(annot_cols, pval_df)
  
    part_batch_data[[batch_idx]] <- final_df
    cat(sprintf("Part %02d - batch %03d finished. Waiting merging.\n", 
                part, batch_idx))

    rm(annot_cols, g_data, K_PRECOMP, precomputed_data, results, test_data, Qk_t, G_rot)
    rm(p_value_output2, combined_p, pval_df, final_df)
    gc()
  }

  # After batch processing is complete, merge all batch data of the current part by row.
  if (length(part_batch_data) > 0) {
    
    part_combined_data <- rbindlist(part_batch_data, use.names = TRUE)
    
    fwrite(part_combined_data, output_file)
    cat(sprintf("Part %02d, %d batches merged. Saved: %s\n", 
                part, total_batches, output_file))
    
    rm(part_combined_data, part_batch_data)
    gc()
  } else {
    warning(sprintf("Part %02d has no information. No output.", part))
  }

  cat(sprintf("Part %02d %d SNP batched finished. \n", part, total_batches))
}


# ---------------------------------------------------------------------------
# 2. Statistics
# ---------------------------------------------------------------------------
library(data.table)
fisher_merge_p <- function(p_vec) {
  if (any(is.na(p_vec))) return(NA_real_)
  -2 * sum(log(p_vec))
}


out_dir <- "./rtMLMM/data/real_data/human"

parts <- sprintf("%02d", 1:33)
annot_cols <- c("SNP", "CHR", "POS")
top2_col <- "Top2" 
block_size <- 65000

suffix <- "quantile"  

input_tpl  <- sprintf("result_larger_%s_part%%s.csv", suffix)#_10
output_file <- file.path(
  out_dir,
  sprintf("SNP_Top2_Fisher_T_Summary_%s.csv", suffix)
)
first_file <- file.path(out_dir, sprintf(input_tpl, parts[1]))
dt_base <- fread(first_file, select = annot_cols)

n_snps   <- nrow(dt_base)
n_blocks <- ceiling(n_snps / block_size)

message("Total SNP rows: ", n_snps)
message("Block size     : ", block_size)
message("Total blocks   : ", n_blocks)

message("Reading all Top2 columns ...")

top2_list <- vector("list", length(parts))
names(top2_list) <- paste0("Top2_part", parts)

for (i in seq_along(parts)) {
  file_path <- file.path(out_dir, sprintf(input_tpl, parts[i]))
  
  top2_list[[i]] <- fread(
    file_path,
    select = top2_col
  )[[1]]
}

top2_matrix_all <- do.call(cbind, top2_list)
colnames(top2_matrix_all) <- paste0("Top2_part", parts)

stopifnot(nrow(top2_matrix_all) == n_snps)

message("Top2 matrix loaded: ",
        nrow(top2_matrix_all), " x ", ncol(top2_matrix_all))


# Omnibus Statistics
for (b in seq_len(n_blocks)) {

  idx <- ((b - 1) * block_size + 1):min(b * block_size, n_snps)
  
  block_base <- dt_base[idx]
  block_top2 <- top2_matrix_all[idx, , drop = FALSE]

  block_res <- cbind(
    block_base,
    as.data.table(block_top2)
  )

  top2_cols <- colnames(block_top2)

  block_res[, T := apply(.SD, 1, fisher_merge_p),
            .SDcols = top2_cols]

  fwrite(block_res, output_file, append = b != 1)

  message("Block ", b, "/", n_blocks,
          " done. Rows: ", nrow(block_res))

  rm(block_res)
  gc()
}

message("All done!")
message("Output file: ", output_file)



# ---------------------------------------------------------------------------
# 3. Find approximately independent SNPs
# ---------------------------------------------------------------------------

library(data.table)

bim_file <- "./rtMLMM/data/real_data/human/geno4_larger.bim"
out_file <- "./rtMLMM/data/real_data/human/independent_snps_final.txt"

bim <- fread(bim_file, header = FALSE)
setnames(bim, c("chr","snp","cm","pos","a1","a2"))

window_size <- 1e5    # 1e5 bp

# Each chromosome
indep_snps_list <- vector("list", length = length(unique(bim$chr)))
chr_ids <- sort(unique(bim$chr))

for(i in seq_along(chr_ids)){
  ch <- chr_ids[i]
  chr_snps <- bim[chr == ch]
  
  last_pos <- chr_snps$pos[1]
  selected_idx <- 1   # Keep the first SNP
  
  for(j in 2:nrow(chr_snps)){
    if(chr_snps$pos[j] - last_pos >= window_size){
      selected_idx <- c(selected_idx, j)
      last_pos <- chr_snps$pos[j]
    }
  }
  
  indep_snps_list[[i]] <- chr_snps[selected_idx, ]

  cat("Chromosome", ch, "independent SNPs:", length(selected_idx), "\n")
}

# Merge all chromosomes
# --------------------------
independent_snps <- rbindlist(indep_snps_list)
fwrite(independent_snps, out_file, col.names = TRUE, quote = FALSE)
cat("Number of approximately independent SNPs:", nrow(independent_snps), "\n")
cat("Approximately independent SNPs Saved:", out_file, "\n")

# ---------------------------------------------------------------------------
# 4. Distribution Approximation via GGD
# ---------------------------------------------------------------------------
library(data.table)
library(moments)
source("./rtMLMM/code/Testmatrix_pcomb.r") 

out_dir <- "./rtMLMM/data/real_data/human/"
norm_type <- "quantile" 
out_file <- "./rtMLMM/data/real_data/human/independent_snps_final.txt"

summary_file <- file.path(out_dir, sprintf("SNP_Top2_Fisher_T_Summary_%s.csv", norm_type))
total_snp <- fread(summary_file)

indep_snps <- fread(out_file, header = TRUE)    # col name: chr snp cm pos a1 a2

subset_snp <- total_snp[SNP %in% indep_snps$snp, ]

cat("Approximately independent SNP count:", nrow(subset_snp), "\n")

T_valid <- subset_snp[!is.na(T), T]

cat("Number of SNPs with non-missing statistics:", length(T_valid), "\n")

lower <- quantile(T_valid, 0.005)
upper <- quantile(T_valid, 0.995)
T_clean <- T_valid[T_valid >= lower & T_valid <= upper]
cat("Valid SNP count:", length(T_clean), "\n")

est <- fit.generalgamma.mle(T_clean, inits = NULL)
mu    <- est["mu"]
sigma <- est["sigma"]
Q     <- est["Q"]
total_snp[, pgengamma := sapply(T, function(stat) {
  if (is.na(stat)) return(9)
  pgengamma(q =stat, mu=mu, sigma=sigma, Q=Q, lower.tail = FALSE)
})]

output_pvalue_file <- file.path(out_dir, sprintf("SNP_Top2_Fisher_T_Summary_with_Pvalue_bp_%s.csv", norm_type))
fwrite(total_snp, output_pvalue_file)
cat("p-value calculation finished.\n")


# ======================================================================================
# step4. Distribution fit plotting and real data analysis plotting
# ======================================================================================


# ---------------------------------------------------------------------------
# 1. Comparison of distribution fitting effects
# ---------------------------------------------------------------------------
library(data.table)
library(ggplot2)
library(flexsurv)
library(fitdistrplus)
library(moments)
library(patchwork)    
source("./rtMLMM/code/Testmatrix_pcomb.r") 
out_dir <- "./rtMLMM/data/real_data/human/"
norm_type <- "quantile"
out_file <- "./rtMLMM/data/real_data/human/independent_snps_final.txt"

summary_file <- file.path(out_dir, sprintf("SNP_Top2_Fisher_T_Summary_%s.csv", norm_type))
total_snp <- fread(summary_file)
indep_snps <- fread(out_file, header = TRUE)
subset_snp <- total_snp[SNP %in% indep_snps$snp, ]

cat("Approximately independent SNP count:", nrow(subset_snp), "\n")
T_valid <- subset_snp[!is.na(T), T]
cat("Number of SNPs with non-missing statistics:", length(T_valid), "\n")

lower <- quantile(T_valid, 0.005)
upper <- quantile(T_valid, 0.995)
T_clean <- T_valid[T_valid >= lower & T_valid <= upper]
cat("Valid SNP count:", length(T_clean), "\n")

# Zoom in on the tail section
tail_quantile <- 0.98
tail_threshold <- quantile(T_clean, tail_quantile)

font_family <- if (.Platform$OS.type == "windows") "Times New Roman" else "Times"

df <- data.frame(T = T_clean)

#  Distribution Approximation
est_gg <- fit.generalgamma.mle(T_clean, inits = NULL)
mu_gg    <- as.numeric(est_gg["mu"])    
sigma_gg <- as.numeric(est_gg["sigma"]) 
Q_gg     <- as.numeric(est_gg["Q"])     
gg_param_label <- paste0("Generalized Gamma (", round(mu_gg,2), ", ", round(sigma_gg,2), ", ", round(Q_gg,2), ")")
cat("=== 1. Generalized Gamma ===\n")
cat("mu:", mu_gg, "\n")
cat("sigma:", sigma_gg, "\n")
cat("Q:", Q_gg, "\n\n")
x_vals <- seq(min(T_clean), max(T_clean), length.out = 1000)
y_gg <- dgengamma(x_vals, mu = mu_gg, sigma = sigma_gg, Q = Q_gg)
df_gg <- data.frame(x = x_vals, y = y_gg, Distribution = "Generalized Gamma")


fit_gamma <- fitdistr(T_clean, "gamma")
shape_gamma <- as.numeric(fit_gamma$estimate["shape"]) 
rate_gamma  <- as.numeric(fit_gamma$estimate["rate"])  
scale_gamma <- 1 / rate_gamma                          
gamma_param_label <- paste0("Gamma (", round(shape_gamma,2), ", ", round(scale_gamma,2), ")")
cat("=== 2. Gamma ===\n")
cat("shape:", shape_gamma, "\n")
cat("rate:", rate_gamma, "\n")
cat("scale:", scale_gamma, "\n\n")
y_gamma <- dgamma(x_vals, shape = shape_gamma, scale = scale_gamma) 
df_gamma <- data.frame(x = x_vals, y = y_gamma, Distribution = "Gamma")


fit_lognorm <- fitdistr(T_clean[T_clean > 0], "lognormal")
mu_lognorm  <- as.numeric(fit_lognorm$estimate["meanlog"]) 
sigma_lognorm <- as.numeric(fit_lognorm$estimate["sdlog"])  
lognorm_param_label <- paste0("Log-normal (", round(mu_lognorm,2), ", ", round(sigma_lognorm,2), ")")
cat("=== 3. Log-normal ===\n")
cat("meanlog:", mu_lognorm, "\n")
cat("sdlog:", sigma_lognorm, "\n\n")
y_lognorm <- dlnorm(x_vals, meanlog = mu_lognorm, sdlog = sigma_lognorm) 
df_lognorm <- data.frame(x = x_vals, y = y_lognorm, Distribution = "Log-normal")


negloglik_chisq <- function(df, x) {
  -sum(dchisq(x, df = df, log = TRUE))
}
m <- mean(T_clean)
opt <- optimize(negloglik_chisq, interval = c(1e-3, 10 * m), x = T_clean)
df_chisq <- as.numeric(opt$minimum) 
chisq_param_label <- paste0("Chi-square (", round(df_chisq,2), ")")
cat("=== 4. Chi-square ===\n")
cat("df:", df_chisq, "\n\n")
y_chisq <- dchisq(x_vals, df = df_chisq) 
df_chisq_plot <- data.frame(x = x_vals, y = y_chisq, Distribution = "Chi-square")


fit_weibull <- fitdistr(T_clean, "weibull")
shape_weibull <- as.numeric(fit_weibull$estimate["shape"])  
scale_weibull <- as.numeric(fit_weibull$estimate["scale"])  
weibull_param_label <- paste0(
  "Weibull (",
  round(shape_weibull, 2), ", ",
  round(scale_weibull, 2), ")"
)

cat("=== 3. Weibull ===\n")
cat("shape:", shape_weibull, "\n")
cat("scale:", scale_weibull, "\n\n")

y_weibull <- dweibull(x_vals, shape = shape_weibull, scale = scale_weibull)

df_weibull <- data.frame(
  x = x_vals,
  y = y_weibull,
  Distribution = "Weibull"
)

library(goftest)
library(philentropy)
library(dplyr)
library(flexsurv)
n <- length(T_clean)
x_vals <- seq(min(T_clean), max(T_clean), length.out = 2000)
x_step <- diff(x_vals)[1] 

dists_list <- list(
  "Generalized Gamma" = list(
    pdf = function(x) dgengamma(x, mu = mu_gg, sigma = sigma_gg, Q = Q_gg),
    cdf = function(x) pgengamma(x, mu = mu_gg, sigma = sigma_gg, Q = Q_gg),
    nparams = 3 
  ),
  "Gamma" = list(
    pdf = function(x) dgamma(x, shape = shape_gamma, scale = scale_gamma),
    cdf = function(x) pgamma(x, shape = shape_gamma, scale = scale_gamma),
    nparams = 2
  ),
  "Log-normal" = list(
    pdf = function(x) dlnorm(x, meanlog = mu_lognorm, sdlog = sigma_lognorm),
    cdf = function(x) plnorm(x, meanlog = mu_lognorm, sdlog = sigma_lognorm),
    nparams = 2
  ),
  "Chi-square" = list(
    pdf = function(x) dchisq(x, df = df_chisq),
    cdf = function(x) pchisq(x, df = df_chisq),
    nparams = 1
  ),
  "Weibull" = list(
    pdf = function(x) dweibull(x, shape = shape_weibull, scale = scale_weibull),
    cdf = function(x) pweibull(x, shape = shape_weibull, scale = scale_weibull),
    nparams = 2
  )
)

# AD distance
ad_stat_manual <- function(x, cdf_fun) {
  x_sorted <- sort(x)
  n <- length(x_sorted)
  F_x <- cdf_fun(x_sorted)
  F_x[F_x == 0] <- 1e-10
  F_x[F_x == 1] <- 1 - 1e-10
  ad <- -n - (1/n) * sum((2*1:n - 1) * (log(F_x) + log(1 - rev(F_x))))
  return(ad)
}

# Fitting metrics
results <- data.frame(
  Distribution = character(),
  AIC = numeric(),
  BIC = numeric(),
  KS_stat = numeric(),
  AD_stat = numeric(),
  KL_divergence = numeric(),
  stringsAsFactors = FALSE
)
dens_emp <- density(T_clean, bw = "SJ") 
p_emp <- approx(dens_emp$x, dens_emp$y, xout = x_vals, rule = 2)$y
p_emp[is.na(p_emp)] <- 0  
p_emp <- p_emp + 1e-10    

p_emp <- p_emp / (sum(p_emp) * x_step)

for (dist_name in names(dists_list)) {
  dist_info <- dists_list[[dist_name]]
  pdf_fun <- dist_info$pdf
  cdf_fun <- dist_info$cdf
  k <- dist_info$nparams  # 
  
  # -----------------
  # 1. AIC/BIC
  # -----------------
  pdf_vals <- pdf_fun(T_clean)
  pdf_vals[pdf_vals <= 0] <- 1e-10   
  loglik <- sum(log(pdf_vals))
  aic <- -2*loglik + 2*k
  bic <- -2*loglik + log(n)*k
  
  # -----------------
  # 2. KS (Deprecated)
  # -----------------
  ks_result <- ks.test(T_clean, cdf_fun)
  ks_stat <- as.numeric(ks_result$statistic)
  
  # -----------------
  # 3. AD (Deprecated)
  # -----------------
  ad_stat <- ad_stat_manual(T_clean, cdf_fun)

  # -----------------
  # 3. KL
  # -----------------  
  p_model <- pdf_fun(x_vals)
  p_model[p_model <= 0] <- 1e-10   
  p_model <- p_model / (sum(p_model) * x_step)
  kl_div <- sum(p_emp * log(p_emp / p_model)) * x_step
  

  results <- rbind(results, data.frame(
    Distribution = dist_name,
    AIC = aic,
    BIC = bic,
    KS_stat = ks_stat,
    AD_stat = ad_stat,
    KL_divergence = kl_div
  ))
}

print(results)

p_main <- ggplot(df, aes(x = T)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 50,
    fill = "lightblue",
    color = "black",
    alpha = 0.6
  ) +


  geom_line(
    data = df_gg,
    aes(x = x, y = y, color = Distribution),
    linewidth = 1.2
  ) +
  geom_line(data = df_gamma,
            aes(x = x, y = y, color = Distribution),
            linewidth = 1.2, linetype = "dotdash") +
  geom_line(data = df_lognorm,
            aes(x = x, y = y, color = Distribution),
            linewidth = 1.2, linetype = "dotdash") +
  geom_line(data = df_chisq_plot,
            aes(x = x, y = y, color = Distribution),
            linewidth = 1.2, linetype = "dotdash") +
  geom_line(data = df_weibull,
            aes(x = x, y = y, color = Distribution),
            linewidth = 1.2, linetype = "dotdash") +

  labs(
    title = " ",
    x = "TOP2-O",
    y = "Density",
    color = "Fitted Distribution"
  ) +

  scale_color_manual(
    values = c(
      "Generalized Gamma" = "red",
      "Gamma"             = "blue",
      "Log-normal"        = "green",
      "Chi-square"        = "purple",
      "Weibull"           = "orange"
    ),
    breaks = c(
      "Generalized Gamma",
      "Gamma",
      "Log-normal",
      "Chi-square",
      "Weibull"
    ),
    labels = c(
      gg_param_label,
      gamma_param_label,
      lognorm_param_label,
      chisq_param_label,
      weibull_param_label
    )
  ) +

  theme(
    text = element_text(family = font_family),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8, fill = NA),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text  = element_text(size = 11, face = "bold", color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "inside",
    legend.position.inside = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text  = element_text(size = 10)
  )
hist_full_data <- ggplot_build(
  ggplot(df, aes(x = T)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50)
)$data[[1]]

hist_tail_data <- hist_full_data[hist_full_data$x >= tail_threshold, ]

gg_tail       <- df_gg[df_gg$x >= tail_threshold, ]
gamma_tail    <- df_gamma[df_gamma$x >= tail_threshold, ]
lognorm_tail  <- df_lognorm[df_lognorm$x >= tail_threshold, ]
chisq_tail    <- df_chisq_plot[df_chisq_plot$x >= tail_threshold, ]
weibull_tail  <- df_weibull[df_weibull$x >= tail_threshold, ]

p_inset <- ggplot() +
  geom_col(
    data = hist_tail_data,
    aes(x = x, y = density),
    fill = "lightblue",
    color = "black",
    alpha = 0.6,
    width = diff(hist_full_data$x)[1]
  ) +

  geom_line(data = gg_tail,
            aes(x = x, y = y, color = Distribution),
            linewidth = 1) +
  geom_line(data = gamma_tail,
            aes(x = x, y = y, color = Distribution),
            linewidth = 1, linetype = "dotdash") +
  geom_line(data = lognorm_tail,
            aes(x = x, y = y, color = Distribution),
            linewidth = 1, linetype = "dotdash") +
  geom_line(data = chisq_tail,
            aes(x = x, y = y, color = Distribution),
            linewidth = 1, linetype = "dotdash") +
  geom_line(data = weibull_tail,
            aes(x = x, y = y, color = Distribution),
            linewidth = 1, linetype = "dotdash") +

  scale_x_continuous(
    limits = c(tail_threshold, max(T_clean)),
    oob = scales::oob_keep
  ) +

  labs(x = "TOP2-O (tail)", y = "Density") +

  scale_color_manual(
    values = c(
      "Generalized Gamma" = "red",
      "Gamma"             = "blue",
      "Log-normal"        = "green",
      "Chi-square"        = "purple",
      "Weibull"           = "orange"
    )
  ) +

  guides(color = "none") +

  theme(
    text = element_text(family = font_family, size = 8),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.6, fill = NA),
    axis.title = element_text(size = 9, face = "bold"),
    axis.text  = element_text(size = 8, face = "bold", color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    plot.margin = margin(1, 1, 1, 1, "mm")
  )
  
p_final <- p_main +
  inset_element(
    p = p_inset,
    left = 0.5,
    bottom = 0.15,
    right = 0.95,
    top = 0.65,
    align_to = "panel"
  )

ggsave(file.path(out_dir, "T_histogram_with_tail_inset_4dist.png"), plot = p_final, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "T_histogram_with_tail_inset_4dist.pdf"), plot = p_final, width = 8, height = 6, device = cairo_pdf)
cat("The final chart has been saved.\n")



# ---------------------------------------------------------------------------
# 2. QQ\Manhattan\Bubble
# ---------------------------------------------------------------------------
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(scales)
library(patchwork)
library(stringr)  

out_dir <- "./rtMLMM/data/real_data/human"
norm_type <- "quantile"
output_dir <- "./rtMLMM/data/real_data/human"


font_family <- "serif"    
axis_title_size <- 14
axis_text_size <- 12
plot_title_size <- 14
legend_text_size <- 12
legend_title_size <- 12

# ==================== 1. QQ ====================
read_file <- file.path(out_dir, sprintf("SNP_Top2_Fisher_T_Summary_with_Pvalue_bp_%s.csv", norm_type))
dt <- fread(read_file)
pvals <- dt$pgengamma
pvals <- pvals[!is.na(pvals) & pvals <= 1 & pvals > 0]

n <- length(pvals)
observed <- -log10(sort(pvals))
i <- 1:n
expected_p <- (i - 0.5) / n
expected <- -log10(expected_p)
qq_data <- data.table(expected = expected, observed = observed)

p_qq <- ggplot(qq_data, aes(x = expected, y = observed)) +
  geom_point(size = 0.8, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "solid") +
  labs(
    x = bquote(Expected - log[10]*textstyle("(")*p*textstyle(")")),
    y = bquote(Observed - log[10]*textstyle("(")*p*textstyle(")")),
    title = "A"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(family = font_family, face = "bold", size = axis_title_size),
    axis.text = element_text(family = font_family, face = "bold", size = axis_text_size),
    plot.title = element_text(family = font_family, face = "bold", size = plot_title_size, hjust = 0),
    panel.grid = element_blank(),
    legend.position = "none"
  )

# ==================== 2. Manhattan ====================
plot_data <- dt[!is.na(pgengamma) & pgengamma > 0 & pgengamma <= 1, ]
plot_data[, logP := -log10(pgengamma)]

chr_lengths <- plot_data[, .(max_pos = max(POS)), by = CHR][order(as.numeric(CHR))]
chr_offsets <- cumsum(c(0, chr_lengths$max_pos[-length(chr_lengths$max_pos)]))
names(chr_offsets) <- chr_lengths$CHR
plot_data[, POS_cum := POS + chr_offsets[CHR]]

chr_centers <- plot_data[, .(center = mean(POS_cum)), by = CHR][order(as.numeric(CHR))]
chr_colors <- rep(c("gray70", "gray40"), length.out = nrow(chr_lengths))
names(chr_colors) <- chr_lengths$CHR
plot_data[, col := chr_colors[CHR]]
sig_threshold <- -log10(5e-8)

p_mh <- ggplot(plot_data, aes(x = POS_cum, y = logP)) +
  geom_point(aes(color = col), alpha = 0.6, size = 1.2) +
  scale_color_identity() +
  scale_x_continuous(breaks = chr_centers$center, labels = chr_centers$CHR, expand = c(0, 0)) +
  scale_y_continuous(breaks = pretty_breaks(n = 4), expand = c(0, 0.5)) +
  geom_hline(yintercept = sig_threshold, linetype = "dashed", color = "red", linewidth = 0.5) +
  labs(x = "Chromosome", y = bquote(-log[10]*textstyle("(")*p*textstyle(")")), title = "B") +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(family = font_family, face = "bold", size = axis_title_size),
    axis.text = element_text(family = font_family, face = "bold", size = axis_text_size),
    plot.title = element_text(family = font_family, face = "bold", size = plot_title_size, hjust = 0),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none"
  )

# ==================== 3. Bubble ====================
bubble_file <- "./rtMLMM/data/real_data/human/WikiPathways_2024_Human_table.txt"
df <- fread(bubble_file)

df_top10 <- df %>%
  separate(Overlap, into = c("gene_count", "gene_background"), sep = "/", convert = TRUE) %>%
  mutate(
    gene_ratio = gene_count / gene_background,
    neg_log10_p = -log10(`P-value`),

    Term_clean = str_remove(Term, "\\s+WP.*$")  
  ) %>%
  arrange(desc(neg_log10_p)) %>%
  slice_head(n = 10) %>%
  mutate(Term_clean = factor(Term_clean, levels = rev(Term_clean)))  

bubble_plot <- ggplot(df_top10, aes(
  x = gene_ratio, y = Term_clean,  
  size = gene_count, color = neg_log10_p
)) +
  geom_point(alpha = 0.9, stroke = 0.5) +
  scale_size_continuous(range = c(5, 20), name = "Gene counts") +
  scale_color_gradientn(colors = c("blue", "purple", "red"), name = "-log10(P Value)") +
  labs(x = "Gene Ratio", y = "", title = "C") +
  theme_bw() +
  theme(
    axis.title = element_text(family = font_family, face = "bold", size = axis_title_size),
    axis.text = element_text(family = font_family, face = "bold", size = axis_text_size),
    plot.title = element_text(family = font_family, face = "bold", size = plot_title_size, hjust = 0),
    legend.text = element_text(family = font_family, face = "bold", size = legend_text_size),
    legend.title = element_text(family = font_family, face = "bold", size = legend_title_size),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

top_row <- p_qq | bubble_plot

final_plot <- top_row / p_mh + 
              plot_layout(heights = c(1, 1)) & 
              theme(text = element_text(family = font_family))

ggsave(file.path(output_dir, "Combined_QQ_Manhattan_Bubble_Clean.pdf"),
       plot = final_plot, width = 16, height = 12, device = "pdf")

ggsave(file.path(output_dir, "Combined_QQ_Manhattan_Bubble_Clean.png"),
       plot = final_plot, width = 16, height = 12, dpi = 300)
       
       
# ======================================================================================
# step5. Summarize all results
# ======================================================================================

# ---------------------------------------------------------------------------
# 1. GWAS
# ---------------------------------------------------------------------------
library(data.table)

OUTPUT_DIR <- "./rtMLMM/data/real_data/human/"
plink_prefix <- paste0(OUTPUT_DIR, "geno4_larger")    # geno4_larger.bed/.bim/.fam


final_df <- fread(
  file = paste0(OUTPUT_DIR, "SNP_Top2_Fisher_T_Summary_with_Pvalue_bp_quantile.csv"),
  stringsAsFactors = FALSE
)

bim_data <- fread(
  file = paste0(plink_prefix, ".bim"),
  col.names = c("CHR", "SNP", "CM", "POS", "A1", "A2"),
  stringsAsFactors = FALSE
)  
bim_data <- bim_data[, .(CHR, POS, SNP, A1, A2)]


system(paste0(
  "plink --bfile ", plink_prefix, 
  " --freq --out ", plink_prefix, "_freq"
))

frq_data <- fread(
  file = paste0(plink_prefix, "_freq.frq"),
  stringsAsFactors = FALSE
)

frq_data <- frq_data[, .(CHR, SNP, A1, A2, MAF)]


snp_info <- merge(
  bim_data[, .(CHR, POS, SNP)],  
  frq_data[, .(CHR, SNP, A1, A2, MAF)], 
  by = c("CHR", "SNP"),
  all.x = TRUE 
)
snp_info[, `:=`(
  major_allele = A2,  
  minor_allele = A1   
)]

snp_info[, c("A1", "A2") := NULL]

final_df[, CHR := as.integer(CHR)]
final_df[, POS := as.integer(POS)]
snp_info[, CHR := as.integer(CHR)]
snp_info[, POS := as.integer(POS)]

result_df <- merge(
  final_df[, .(CHR, POS, SNP, pgengamma)],  
  snp_info[, .(CHR, POS, SNP, major_allele, minor_allele, MAF)],
  by = c("CHR", "POS", "SNP"),
  all.x = TRUE  
)

result_df <- result_df[, .(
  CHR,          # chromosome
  POS,          # position
  SNP,          # SNP name
  major_allele,
  minor_allele,
  MAF,
  pgengamma     # p-value
)]
cat("===== Missing values =====\n")

major_na_before <- sum(is.na(result_df$major_allele))
cat("Number of missing major_allele:", major_na_before, "\n")

minor_na_before <- sum(is.na(result_df$minor_allele))
cat("Number of missing minor_allele:", minor_na_before, "\n")

maf_na_before <- sum(is.na(result_df$MAF))
cat("Number of missing MAF:", maf_na_before, "\n")
result_df[is.na(major_allele), major_allele := "NA"]
result_df[is.na(minor_allele), minor_allele := "NA"]
result_df[is.na(MAF), MAF := 0]

cat("\n===== Fill missing values =====\n")

major_na_after <- sum(result_df$major_allele == "NA")
cat("Number of major_alleles filled with NA:", major_na_after, "\n")
minor_na_after <- sum(result_df$minor_allele == "NA")
cat("Number of minor_allele filled with NA", minor_na_after, "\n")
maf_zero_after <- sum(result_df$MAF == 0)
cat("Number of MAF filled with zero", maf_zero_after, "\n")

total_rows <- nrow(result_df)
cat("\n Number of missing rows", total_rows, "\n")
cat("Missing rate:", round((major_na_before / total_rows) * 100, 2), "%\n")

output_file <- paste0(OUTPUT_DIR, "SNP_with_allele_MAF_pgengamma_final.csv")
fwrite(result_df, file = output_file, row.names = FALSE, na = "")
cat("\n Final file saved.", output_file, "\n")
cat("Preview first 5 rows:\n")
print(result_df[1:5, ])


# ---------------------------------------------------------------------------
# 2. Significant SNP
# ---------------------------------------------------------------------------

library(data.table)

OUTPUT_DIR <- "./rtMLMM/data/real_data/human/"
input_file <- paste0(OUTPUT_DIR, "SNP_with_allele_MAF_pgengamma_final.csv")
output_filtered_file <- paste0(OUTPUT_DIR, "SNP_with_allele_MAF_pgengamma_significant.csv")

if (!file.exists(input_file)) {
  stop("Error. ", input_file, " not exist")
}

df <- fread(input_file, stringsAsFactors = FALSE)

cat("===== colname =====\n")
print(colnames(df))
cat("\n===== preview first 3 rows =====\n")
print(df[1:3, ])

df[, pgengamma := as.numeric(pgengamma)]

filtered_df <- df[!is.na(pgengamma) & pgengamma < 5e-8, ]

total_rows <- nrow(df)
filtered_rows <- nrow(filtered_df)

cat("\n===== Screening for significant SNPs =====\n")
cat("Number of raw snp:", total_rows, "\n")
cat("Number of significant snp (pgengamma < 5e-8):", filtered_rows, "\n")
cat("Rate:", round((filtered_rows / total_rows) * 100, 4), "%\n")

fwrite(filtered_df, file = output_filtered_file, row.names = FALSE, na = "")

cat("Significant snp information saved:", output_filtered_file, "\n")

if (filtered_rows > 0) {
  cat("\n Preview first 5 significant snp rows:\n")
  print(filtered_df[1:5, ])
} else {
  cat("\n no pgengamma < 5e-8 snp\n")
}

# ---------------------------------------------------------------------------
# 3. Significant gene (SNP mapping to genes)
# ---------------------------------------------------------------------------
# [Shell Command] The following steps were performed in Linux terminal:
# use mgmma
# wget -P ./rtMLMM/data/real_data/human/ https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.gene.loc
# ./magma \
#   --annotate window=5,5 \
#   --snp-loc ./rtMLMM/data/real_data/human/geno4_larger.bim \
#   --gene-loc ./rtMLMM/data/real_data/human/NCBI37.3.gene.loc \
#   --out ./rtMLMM/data/real_data/human/magma_annot_larger
# ---------------------------------------------------------------------------

library(data.table)
library(dplyr)

OUTPUT_DIR <- "./rtMLMM/data/real_data/human/"

# ============================
# 1) significant SNP
# ============================
snp_sig <- fread(file.path(
  OUTPUT_DIR,
  "SNP_with_allele_MAF_pgengamma_significant.csv"
))

lines <- readLines(file.path(
  OUTPUT_DIR,
  "magma_annot_larger.genes.annot"
))
lines <- lines[!grepl("^#", lines)]


magma_list <- lapply(lines, function(line) {
  parts <- strsplit(line, "\\s+")[[1]]
  if (length(parts) < 3) return(NULL)

  gene_id <- parts[1]
  chr_coords <- parts[2]  # CHR:start:end
  chr_split <- strsplit(chr_coords, ":")[[1]]

  CHR <- chr_split[1]
  gene_start <- as.numeric(chr_split[2])
  gene_end   <- as.numeric(chr_split[3])

  snps <- parts[-c(1, 2)]
  if (length(snps) == 0) return(NULL)

  data.frame(
    gene_id    = gene_id,
    CHR        = CHR,
    gene_start = gene_start,
    gene_end   = gene_end,
    SNP        = snps,
    stringsAsFactors = FALSE
  )
})

magma_long <- bind_rows(magma_list) %>%
  filter(!is.na(SNP), SNP != "")

magma_sig <- magma_long %>%
  filter(SNP %in% snp_sig$SNP)


snp_pos <- snp_sig %>% select(SNP, BP = POS)
magma_sig <- magma_sig %>%
  left_join(snp_pos, by = "SNP")


magma_sig <- magma_sig %>%
  mutate(
    dist = ifelse(
      BP < gene_start, gene_start - BP,
      ifelse(BP > gene_end, BP - gene_end, 0)
    )
  )

# At most one gene is retained for each SNP.
snp_to_gene <- magma_sig %>%
  group_by(SNP) %>%
  slice_min(order_by = dist, n = 1, with_ties = FALSE) %>%
  ungroup()

gene_summary <- snp_to_gene %>%
  group_by(CHR, gene_start, gene_end, gene_id) %>%
  summarise(
    n_SNP = n(),
    SNP_list = paste(unique(SNP), collapse = ","),
    .groups = "drop"
  )


gene_loc <- fread(file.path(
  OUTPUT_DIR,
  "NCBI37.3.gene.loc"
))
colnames(gene_loc) <- c(
  "gene_id", "CHR", "gene_start", "gene_end", "strand", "gene_name"
)

gene_loc <- gene_loc %>%
  mutate(gene_id = as.character(gene_id))

gene_summary <- gene_summary %>%
  mutate(gene_id = as.character(gene_id)) %>%
  left_join(
    gene_loc %>% select(gene_id, official_gene_name = gene_name),
    by = "gene_id"
  ) %>%
  mutate(gene_name = official_gene_name) %>%
  select(CHR, gene_start, gene_end, gene_name, n_SNP, SNP_list)

chr_levels <- c(as.character(1:22))

gene_summary <- gene_summary %>%
  mutate(
    CHR = factor(CHR, levels = chr_levels)
  ) %>%
  arrange(CHR, gene_start)

fwrite(
  gene_summary,
  file.path(
    OUTPUT_DIR,
    "gene_significant_summary_gene_centric.csv"
  )
)

cat(
  "Done! Results saved: gene_significant_summary_gene_centric.csv\n"
)

# ---------------------------------------------------------------------------
# 4. PubMed search
# ---------------------------------------------------------------------------

library(rentrez)
library(data.table)
library(dplyr)
library(xml2)

Entrez_email <- "2684283931@qq.com"
options(timeout = 60)
set.seed(1)

OUTPUT_DIR <- "./rtMLMM/data/real_data/human/"

# phenotype
phenotype_map <- data.frame(
  short = c("ALT","ALB","ALP","ApoA","ApoB","AST","CRP","Ca","TC","Cr","CysC",
            "GGT","Glu","HbA1c","HDL-C","IGF-1","LDL-C","P","SHBG","T","TBil",
            "TP","TG","UA","Urea","Vit D"),
  full = c("Alanine aminotransferase","Albumin","Alkaline phosphatase",
           "Apolipoprotein A","Apolipoprotein B",
           "Aspartate aminotransferase","C-reactive protein","Calcium",
           "Cholesterol","Creatinine","Cystatin C",
           "Gamma glutamyltransferase","Glucose",
           "Glycated haemoglobin","HDL cholesterol","Insulin-like Growth Factor 1",
           "Low Density Lipoprotein Cholesterol","Phosphate","Sex Hormone Binding Globulin",
           "Testosterone","Total bilirubin","Total protein","Triglycerides",
           "Uric Acid","Urea","Vitamin D"),
  stringsAsFactors = FALSE
)

gene_summary <- fread(file.path(
  OUTPUT_DIR,
  "gene_significant_summary_gene_centric.csv"
))
genes <- unique(gene_summary[[4]])
cat("Total genes:", length(genes), "\n")

search_pmids <- function(gene, phenotype_full, retmax = 20) {
  # search "humans[MeSH Terms]"
  query <- paste0(
    gene, "[Title/Abstract] AND ",
    "\"", phenotype_full, "\"[Title/Abstract] AND humans[MeSH Terms]"
  )
  res <- tryCatch(
    entrez_search(
      db = "pubmed",
      term = query,
      retmax = retmax,
      email = Entrez_email
    ),
    error = function(e) NULL
  )
  if (is.null(res) || length(res$ids) == 0) return(character(0))
  unique(res$ids)
}

pmid_to_citation_xml <- function(pmid) {
  xml_txt <- tryCatch(
    entrez_fetch(
      db = "pubmed",
      id = pmid,
      rettype = "xml",
      parsed = FALSE,
      email = Entrez_email
    ),
    error = function(e) NULL
  )
  if (is.null(xml_txt)) return(NA_character_)

  doc <- read_xml(xml_txt)

  author <- xml_find_first(doc, ".//AuthorList/Author[1]/LastName")
  author <- ifelse(length(author) > 0, xml_text(author), NA_character_)

  year <- xml_find_first(doc, ".//PubDate/Year")
  year <- ifelse(length(year) > 0, xml_text(year), NA_character_)

  if (!is.na(author) && !is.na(year)) {
    paste0(author, " et al. (", year, ")")
  } else {
    NA_character_
  }
}

res_list <- list()
idx <- 1

for (g in genes) {
  for (i in seq_len(nrow(phenotype_map))) {

    pheno_short <- phenotype_map$short[i]
    pheno_full  <- phenotype_map$full[i]

    pmids <- search_pmids(g, pheno_full)

    if (length(pmids) > 0) {
      citations <- sapply(pmids, pmid_to_citation_xml)
      res_list[[idx]] <- data.frame(
        Gene = g,
        Phenotype = pheno_short,
        PMID = pmids,
        Citation = citations,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1

      cat(sprintf(
        "Processed: %-8s %-6s PMID=%d\n",
        g, pheno_short, length(pmids)
      ))
    } else {
      cat(sprintf(
        "Processed: %-8s %-6s PMID=0\n",
        g, pheno_short
      ))
    }

    Sys.sleep(0.3) 
  }
}

if(length(res_list) > 0){
  final_res <- bind_rows(res_list)
} else {
  final_res <- data.frame(Gene=character(),
                          Phenotype=character(),
                          PMID=character(),
                          Citation=character(),
                          stringsAsFactors = FALSE)
}

fwrite(
  final_res,
  file.path(OUTPUT_DIR, "gene_phenotype_pubmed_clean_human.csv")
)

cat("DONE: gene_phenotype_pubmed_clean_human.csv\n")

library(data.table)
library(dplyr)
OUTPUT_DIR <- "./rtMLMM/data/real_data/human/"
gene_summary <- fread(
  file.path(OUTPUT_DIR, "gene_significant_summary_gene_centric.csv")
)

pubmed_info <- fread(
  file.path(OUTPUT_DIR, "gene_phenotype_pubmed_clean_human.csv")
)
gene_summary2 <- gene_summary %>%
  rename(Gene = gene_name)

pubmed_first <- pubmed_info %>%
  filter(!is.na(PMID)) %>%      
  group_by(Gene) %>%
  summarise(
    PMID = first(PMID),
    Citation = first(Citation),
    .groups = "drop"
  )

gene_summary_updated <- gene_summary2 %>%
  left_join(pubmed_first, by = "Gene")


out_file <- file.path(
  OUTPUT_DIR,
  "gene_significant_summary_with_pubmed.csv"
)

fwrite(gene_summary_updated, out_file)

cat("====================================\n")
cat("Total genes: ", nrow(gene_summary_updated), "\n")
cat(
  "Genes with PMID: ",
  sum(!is.na(gene_summary_updated$PMID)),
  "\n"
)
cat("Output file: ", out_file, "\n")
cat("DONE.\n")
cat("====================================\n")

