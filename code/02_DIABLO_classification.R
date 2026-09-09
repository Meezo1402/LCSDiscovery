# ==============================================================================
# DIABLO Nested Cross-Validation - HMM Combined-State Classification
# ==============================================================================
# Purpose : Estimate out-of-sample balanced error rate (BER) for a DIABLO
#           (mixOmics::block.splsda) model predicting Phase 3 HMM combined
#           states (Moderate-Severe Affected vs. Unaffected-Mild Affected)
#           from Phase 1 omics blocks (metabolomics, methylation, clinical).
#
# Design  : 5-fold CV repeated 2x (caret::createMultiFolds) = 10 outer test
#           sets. ncomp (via perf) and keepX (via tune.block.splsda) are tuned
#           WITHIN each training fold only, so feature selection and
#           hyperparameter tuning never see the held-out test subjects.
#
# Inputs  : metabo_HMM_ph1static, methyl_HMM_ph1static, clin_HMM_ph1static
#             - subjects x features; rownames are subject IDs
#           key_ph3_combostates / key_ph3static_combostates /
#           key_ph3static_caseVctrl
#             - ID-to-state key dataframes (columns: ID, State)
#           metabo_drop
#             - character vector: metabolomics metadata columns + salicylate
#
# Outputs : BER_values, BER_values_class1, BER_values_class2 (length 10)
#
# Session : R 4.4.1 | mixOmics 6.28.0 | caret 7.0-1 | dplyr 1.1.4 | plyr 1.8.9
#
# Notes   : Metadata columns are dropped positionally from the methylation block
#           (cols 1-3) and clinical block (cols 1, 2, 7, 27, 28, 29).
#           Subject/label alignment relies on all input frames sharing the
#           same subject row order.
# ==============================================================================

library(caret)      # createMultiFolds()
library(mixOmics)   # block.splsda(), perf(), tune.block.splsda()

# ==============================================================================
# 1. Subject IDs and labels
# ==============================================================================
# Subject 99 is excluded.
# NOTE: the exclusion index is computed on key_ph3static_combostates while the
# ID/State vectors come from key_ph3_combostates; this is valid only if the
# two key frames share identical row order. Verify before reuse.

set.seed(123)   # reproducible folds every run

all_ids   <- key_ph3_combostates$ID[-which(key_ph3static_combostates$ID == as.numeric("99"))]
all_label <- key_ph3_combostates$State[-which(key_ph3static_combostates$ID == as.numeric("99"))]

# ==============================================================================
# 2. Cross-validation folds (5-fold x 2 repeats = 10 outer splits)
# ==============================================================================
# createMultiFolds() returns TRAINING indices (not test indices).

set.seed(123)
fold_idx <- caret::createMultiFolds(all_label, k = 5, times = 2)

folds <- lapply(fold_idx, function(train_idx) {
  train_ids <- all_ids[train_idx]
  test_ids  <- setdiff(all_ids, train_ids)
  list(train = train_ids, test = test_ids)
})

# ==============================================================================
# 3. Storage
# ==============================================================================

num_iterations    <- length(folds)   # 10
BER_values        <- numeric(num_iterations)
BER_values_class1 <- numeric(num_iterations)
BER_values_class2 <- numeric(num_iterations)

# ==============================================================================
# 4. Outer CV loop: build blocks, tune within training data, evaluate on test
# ==============================================================================

set.seed(123)
start_time <- Sys.time()

for (i in seq_len(num_iterations)) {
  cat("Fold", i, "\n")

  train_subjects <- folds[[i]]$train
  test_subjects  <- folds[[i]]$test

  # ---- Training blocks (salicylate dropped via metabo_drop; source frames
  #      unchanged) -------------------------------------------------------
  train_data <- list(
    metabo = as.matrix(metabo_HMM_ph1static[rownames(metabo_HMM_ph1static) %in% train_subjects,-c(1:3)]),
    meth   = as.matrix(methyl_HMM_ph1static[rownames(methyl_HMM_ph1static) %in% train_subjects, -c(1:3)]),
    clin   = as.matrix(clin_HMM_ph1static[rownames(clin_HMM_ph1static) %in% train_subjects, -c(1, 2, 7, 27, 28, 29)])
  )

  # ---- Test blocks (same column exclusions as training) -----------------
  test_data <- list(
    metabo = as.matrix(metabo_HMM_ph1static[rownames(metabo_HMM_ph1static) %in% test_subjects, -c(1:3) ]),
    meth   = as.matrix(methyl_HMM_ph1static[rownames(methyl_HMM_ph1static) %in% test_subjects, -c(1:3)]),
    clin   = as.matrix(clin_HMM_ph1static[rownames(clin_HMM_ph1static) %in% test_subjects, -c(1, 2, 7, 27, 28, 29)])
  )

  # ---- Outcome labels (alignment relies on shared subject row order) ----
  Y_train <- key_ph3static_caseVctrl$State[key_ph3static_caseVctrl$ID %in% train_subjects]
  Y_test  <- key_ph3static_caseVctrl$State[key_ph3static_caseVctrl$ID %in% test_subjects]

  # ---- DIABLO design matrix (0.1 off-diagonal balances discrimination
  #      vs. cross-block correlation) -------------------------------------
  des <- matrix(0.1, nrow = length(train_data), ncol = length(train_data),
                dimnames = list(names(train_data), names(train_data)))
  des[row(des) == col(des)] <- 0

  # ---- Tune ncomp within the training fold ------------------------------
  base_model <- block.splsda(train_data, Y_train, ncomp = 5, design = des)
  perf_out   <- perf(base_model, validation = "Mfold", folds = 10, nrepeat = 10)
  ncomp      <- perf_out$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]

  # ---- Tune keepX within the training fold ------------------------------
  keep_grid <- list(metabo = c(50, 100, 200, 300),
                    meth   = c(50, 100, 200, 500),
                    clin   = c(5, 10, 20, 30))

  tune_out <- tune.block.splsda(train_data, Y_train, ncomp = ncomp,
                                test.keepX = keep_grid, design = des,
                                validation = "Mfold", folds = 10, nrepeat = 1,
                                dist = "centroids.dist")

  # ---- Final fold model and held-out evaluation -------------------------
  final_mod <- block.splsda(train_data, Y_train, ncomp = ncomp,
                            keepX = tune_out$choice.keepX, design = des)

  pred <- predict(final_mod, newdata = test_data)
  cm   <- get.confusion_matrix(truth = Y_test,
                               predicted = pred$WeightedVote$centroids.dist[, ncol(pred$WeightedVote$centroids.dist)])

  # Per-class error: rows of cm follow the truth factor-level order;
  # confirm with rownames(cm) on first run.
  BER_values[i]        <- get.BER(cm)
  BER_values_class1[i] <- (1 - diag(cm) / rowSums(cm))[1]
  BER_values_class2[i] <- (1 - diag(cm) / rowSums(cm))[2]
}

end_time <- Sys.time()
tot_time <- round(end_time - start_time, 2)
print(tot_time)

# ==============================================================================
# 5. Summary
# ==============================================================================

cat("\nBER across", num_iterations, "outer folds:\n")
cat("  Overall BER : mean =", round(mean(BER_values), 3),
    "| sd =", round(sd(BER_values), 3), "\n")
cat("  Class 1 err : mean =", round(mean(BER_values_class1), 3),
    "| sd =", round(sd(BER_values_class1), 3), "\n")
cat("  Class 2 err : mean =", round(mean(BER_values_class2), 3),
    "| sd =", round(sd(BER_values_class2), 3), "\n")
