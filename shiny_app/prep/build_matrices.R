################################################################################
# build_matrices.R
#
# Parse the estimated parameter matrices A / R / Lambda / Sigma for every model
# and criterion (pAIC / pBIC) from output/real/comparison_table_<crit>.csv into
# a long, app-friendly table (C-1), plus a model-fit summary and the optimal
# gamma table for context.
#
# The point of C-1: conventional individual models force the cross-variable
# structure to be diagonal (off-diagonal Lambda/Sigma = 0 or absent), while the
# new methodology estimates off-diagonal coupling. The `present` flag
# distinguishes "freely estimated" from "structurally constrained to 0 / not in
# model".
#
# Significance: the table has two SE sources per cell. Psi-based SE is used when
# it exists; where it does not (R and Lambda of the new models) the
# Hessian-based SE is used instead. `signif_source` records which one
# ("psi" / "hessian"; NA when neither SE exists, e.g. Sigma).
#
# App reads only: this script does no estimation.
#
# Output:
#   shiny_app/data/coef_matrices.rds = list(coef = <df>, fit = <df>, gamma = <df>)
#
# Run from the repository root:
#   Rscript shiny_app/prep/build_matrices.R
################################################################################

suppressMessages({
  library(dplyr)
})

REPO_ROOT <- Sys.getenv("MSNR_REPO_ROOT", getwd())
source(file.path(REPO_ROOT, "shiny_app", "R", "constants.R"))

PATH_CT    <- function(crit) {
  file.path(REPO_ROOT, "output", "real", sprintf("comparison_table_%s.csv", crit))
}
PATH_SEL   <- file.path(REPO_ROOT, "output", "real", "model_selection_summary.csv")
PATH_GAMMA <- file.path(REPO_ROOT, "output", "real", "gamma_optimal_summary.csv")
OUT_DIR    <- file.path(REPO_ROOT, "shiny_app", "data")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

cat("=== build_matrices.R (A/R/Lambda/Sigma; app reads only) ===\n")
stopifnot(all(file.exists(vapply(CRITERIA, PATH_CT, character(1)))),
          file.exists(PATH_SEL), file.exists(PATH_GAMMA))

################################################################################
# §1. Parse matrix parameters from comparison_table_<crit>.csv
################################################################################

# NULL / empty / all-NA -> fallback (base R >= 4.4 `%||%` handles NULL only).
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

# Cell text of a significance column: NA -> "".
sig_text <- function(ct, col, i) {
  if (!col %in% colnames(ct)) return("")
  s <- as.character(ct[[col]][i])
  if (is.na(s)) "" else trimws(s)
}
num_at <- function(ct, col, i) {
  if (!col %in% colnames(ct)) return(NA_real_)
  suppressWarnings(as.numeric(ct[[col]][i]))
}

parse_ct <- function(crit) {
  ct <- read.csv(PATH_CT(crit), check.names = FALSE, stringsAsFactors = FALSE)
  param <- ct[[1]]
  # Rows like "Lambda[2,1]"; Sigma is symmetric and stores only [1,1],[1,2],[2,2].
  mrx <- regmatches(param, regexec("^(A|R|Lambda|Sigma)\\[([12]),([12])\\]$", param))
  is_mat <- vapply(mrx, length, integer(1)) == 4
  stopifnot(sum(is_mat) == 15)

  rows <- list()
  for (i in which(is_mat)) {
    mat <- mrx[[i]][2]; r <- as.integer(mrx[[i]][3]); cc <- as.integer(mrx[[i]][4])
    for (mid in MODEL_IDS) {
      est    <- num_at(ct, paste0("estimate.", mid), i)
      se_psi <- num_at(ct, paste0("se_psi.", mid), i)
      se_hes <- num_at(ct, paste0("se_hessian.", mid), i)
      # Choose the SE source by availability of the SE itself (not by whether
      # the significance mark is blank: a blank mark with a Psi-SE means
      # "not significant", and must not fall back to the Hessian).
      src <- if (is.finite(se_psi)) "psi" else if (is.finite(se_hes)) "hessian" else NA_character_
      sig <- switch(src %||% "none",
                    psi     = sig_text(ct, paste0("signif_psi.", mid), i),
                    hessian = sig_text(ct, paste0("signif_hessian.", mid), i),
                    "")

      off_diag <- (r != cc)
      # present = freely estimated. NA -> not in model; off-diagonal exact 0 ->
      # structurally constrained (diagonal restriction of conventional models).
      present <- is.finite(est) && !(off_diag && est == 0)

      rows[[length(rows) + 1L]] <- data.frame(
        model = mid, model_type = unname(MODEL_TYPE[mid]),
        family = unname(MODEL_FAMILY[mid]), criterion = crit,
        matrix = mat, row = r, col = cc,
        estimate = est,
        se = switch(src %||% "none", psi = se_psi, hessian = se_hes, NA_real_),
        signif = sig, signif_source = src, present = present,
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(rows)
}

coef <- bind_rows(lapply(CRITERIA, parse_ct))

# Sigma is symmetric: mirror [1,2] into [2,1] so the 2x2 grid renders fully.
sig12 <- coef %>% filter(matrix == "Sigma", row == 1, col == 2)
stopifnot(nrow(sig12) == length(MODEL_IDS) * length(CRITERIA))
sig21 <- sig12; sig21$row <- 2L; sig21$col <- 1L
coef <- bind_rows(coef, sig21)

################################################################################
# §2. Model-fit summary (context for C-1)
################################################################################

sel <- read.csv(PATH_SEL, stringsAsFactors = FALSE,
                colClasses = c(model_id = "character"))  # keep leading 0 ("0111")
sel$model_id <- as.character(sel$model_id)
stopifnot(setequal(sel$model_id, MODEL_IDS))
fit <- sel %>%
  transmute(
    model = model_id,
    model_type = unname(MODEL_TYPE[model]),
    family = unname(MODEL_FAMILY[model]),
    n_params, d_eff_pAIC, d_eff_pBIC, loglik_pAIC, loglik_pBIC,
    AIC, pAIC, BIC, pBIC, gamma_pAIC, gamma_pBIC,
    R2 = R2_mean, adj_R2 = adj_R2_mean,
    delta_AIC, delta_pAIC, delta_BIC, delta_pBIC
  ) %>%
  arrange(pAIC)

gamma <- read.csv(PATH_GAMMA, stringsAsFactors = FALSE,
                  colClasses = c(model_id = "character"))
gamma$model_id <- as.character(gamma$model_id)
stopifnot(all(gamma$model_id %in% MODEL_IDS))
gamma <- gamma %>% rename(model = model_id) %>%
  mutate(model_type = unname(MODEL_TYPE[model]), .after = model)

saveRDS(list(coef = coef, fit = fit, gamma = gamma),
        file.path(OUT_DIR, "coef_matrices.rds"))

################################################################################
# §3. Sanity checks
################################################################################

cat("=== Sanity checks ===\n")
n_expect <- (15 + 1) * length(MODEL_IDS) * length(CRITERIA)
cat(sprintf("  coef rows: %d (expect %d) | models: %d | criteria: %s | matrices: %s\n",
            nrow(coef), n_expect, n_distinct(coef$model),
            paste(sort(unique(coef$criterion)), collapse = "/"),
            paste(sort(unique(coef$matrix)), collapse = ", ")))
stopifnot(nrow(coef) == n_expect)

cat("  significance source by matrix (present cells):\n")
print(coef %>% filter(present) %>%
        count(matrix, signif_source) %>% as.data.frame())

cat("  off-diagonal Lambda 'present' by family (pAIC):\n")
print(coef %>% filter(criterion == "pAIC", matrix == "Lambda", row != col) %>%
        group_by(family) %>%
        summarise(present = sum(present), total = n(), .groups = "drop") %>%
        as.data.frame())
cat("  A 'present' by model (pAIC; expect 0 for 1101/0101/1001/000d):\n")
a_pres <- coef %>% filter(criterion == "pAIC", matrix == "A") %>%
  group_by(model) %>% summarise(present = sum(present), .groups = "drop")
print(as.data.frame(a_pres[match(MODEL_IDS, a_pres$model), ]))
stopifnot(all(a_pres$present[a_pres$model %in% c("1101", "0101", "1001", "000d")] == 0))

# pAIC vs pBIC estimates: same where gamma is the same in both criteria
d_est <- inner_join(
  coef %>% filter(criterion == "pAIC") %>% select(model, matrix, row, col, a = estimate),
  coef %>% filter(criterion == "pBIC") %>% select(model, matrix, row, col, b = estimate),
  by = c("model", "matrix", "row", "col")) %>%
  group_by(model) %>%
  summarise(max_abs_diff = max(abs(a - b), na.rm = TRUE), .groups = "drop")
cat("  pAIC vs pBIC matrix estimates, max |diff|:\n")
print(as.data.frame(d_est[match(MODEL_IDS, d_est$model), ]))

cat(sprintf("  fit rows: %d (best pAIC: %s %s; best pBIC: %s %s)\n",
            nrow(fit), fit$model[1], fit$model_type[1],
            fit$model[which.min(fit$pBIC)], fit$model_type[which.min(fit$pBIC)]))
cat(sprintf("  gamma rows: %d (%s)\n", nrow(gamma), paste(gamma$model, collapse = ", ")))

cat("\nDone. -> shiny_app/data/coef_matrices.rds\n")
