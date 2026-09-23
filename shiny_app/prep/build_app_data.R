################################################################################
# build_app_data.R
#
# Regenerate the Shiny app's long-format data:
#   - shiny_app/data/app_data.rds   (response / covariate / residual / fitted)
#
# PENALIZED estimation (pAIC / pBIC):
#   Each penalized model has its own optimal gamma per criterion, so fitted
#   values and residuals can differ between pAIC and pBIC. The `criterion`
#   column (pAIC / pBIC; NA for response and covariate rows) keeps both sets.
#
# Residuals are the INNOVATION residuals `y*_resid` (the residuals of the
# likelihood / Sigma estimation). The reduced-form residuals
# `y*_resid_reduced` (= obs - pred) are deliberately NOT read.
# Fitted values are `y*_pred` (reduced-form trend prediction).
#
# Scope:
#   Geography and municipality names do NOT depend on the estimation, so
#   muni_sf.rds and pref_borders.rds are REUSED as-is and are NOT regenerated.
#
# Residual differences (e_A - e_B) are NOT pre-computed: the app computes the
# difference for any model pair (within one criterion) on the fly.
#
# Run from the repository root:
#   Rscript shiny_app/prep/build_app_data.R
################################################################################

suppressMessages({
  library(dplyr)
})

################################################################################
# §0. Paths and parameters
################################################################################

REPO_ROOT <- Sys.getenv("MSNR_REPO_ROOT", getwd())
source(file.path(REPO_ROOT, "shiny_app", "R", "constants.R"))

PATH_DATA    <- file.path(REPO_ROOT, "data", "real", "transformed_data198_for_R.csv")
PATH_MAPPING <- file.path(REPO_ROOT, "data", "real", "region_id_mapping198.csv")
RESID_DIR    <- function(crit) {
  file.path(REPO_ROOT, "output", "real", paste0("fitted_residuals_", crit))
}

OUT_DIR <- file.path(REPO_ROOT, "shiny_app", "data")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

cat("=== build_app_data.R (app_data.rds only; penalized pAIC/pBIC) ===\n")
cat("  REPO_ROOT:", REPO_ROOT, "\n")
stopifnot(file.exists(PATH_DATA), file.exists(PATH_MAPPING),
          all(dir.exists(vapply(CRITERIA, RESID_DIR, character(1)))))

################################################################################
# §1. Mapping: region -> code_5
################################################################################

# Columns: region, code (6-digit), name (Japanese).
mapping <- read.csv(PATH_MAPPING, fileEncoding = "UTF-8")
stopifnot(identical(names(mapping), c("region", "code", "name")))
mapping$code   <- as.character(mapping$code)
mapping$code_5 <- substr(mapping$code, 1, 5)
stopifnot(nrow(mapping) == 198)
cat(sprintf("  Mapping: %d municipalities\n", nrow(mapping)))

# region_id -> code_5 lookup (transformed/residual CSVs key by region)
reg2code <- setNames(mapping$code_5, mapping$region)

################################################################################
# §2. app_data: response + covariate (both years)
################################################################################

cat("=== Building app_data (long) ===\n")
dat <- read.csv(PATH_DATA)
stopifnot(nrow(dat) == 396, setequal(dat$time, as.integer(names(YEARS))))
dat$muni_code <- unname(reg2code[as.character(dat$region)])
dat$year      <- unname(YEARS[as.character(dat$time)])

# variable key -> (column, category, outcome, label), from constants.R
VAR_SPEC <- rbind(
  data.frame(variable = OUTCOMES$key, col = OUTCOMES$outcome,
             category = "response", outcome = OUTCOMES$outcome,
             label = OUTCOMES$label, stringsAsFactors = FALSE),
  data.frame(variable = COVARIATES$key, col = COVARIATES$col,
             category = "covariate", outcome = NA_character_,
             label = COVARIATES$label, stringsAsFactors = FALSE)
)

base_long <- lapply(seq_len(nrow(VAR_SPEC)), function(i) {
  s <- VAR_SPEC[i, ]
  data.frame(
    muni_code = dat$muni_code,
    year      = dat$year,
    category  = s$category,
    variable  = s$variable,
    outcome   = s$outcome,
    model     = NA_character_,
    criterion = NA_character_,
    value     = dat[[s$col]],
    label     = s$label,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

################################################################################
# §3. app_data: residual / fitted for all models x criteria (YEAR_RESID only)
################################################################################

read_resid <- function(crit, mid) {
  fpath <- file.path(RESID_DIR(crit), sprintf("fitted_residuals_%s.csv", mid))
  if (!file.exists(fpath)) stop(sprintf("Missing: %s", fpath))
  df <- read.csv(fpath)
  df <- df[df$time == 2, , drop = FALSE]
  stopifnot(nrow(df) == 198)
  df$muni_code <- unname(reg2code[as.character(df$region)])
  df
}

# Labels carry no criterion: the variable key is shared by pAIC/pBIC and the
# criterion is a separate column (so variable choices do not duplicate).
resid_long <- list()
for (crit in CRITERIA) {
  for (mid in MODEL_IDS) {
    df <- read_resid(crit, mid)
    for (oc in OUTCOMES$outcome) {
      key <- sprintf("%s_%s", OUTCOME_KEY[[oc]], mid)
      for (q in c("residual", "fitted")) {
        col <- paste0(oc, if (q == "residual") "_resid" else "_pred")
        resid_long[[paste(crit, mid, oc, q, sep = "_")]] <- data.frame(
          muni_code = df$muni_code, year = YEAR_RESID, category = q,
          variable = key, outcome = oc, model = mid, criterion = crit,
          value = df[[col]],
          label = sprintf("%s: %s — %s",
                          if (q == "residual") "Residual" else "Fitted",
                          OUTCOME_SHORT[[oc]], MODEL_LABELS[[mid]]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

app_data <- bind_rows(base_long, bind_rows(resid_long))
stopifnot(!anyNA(app_data$muni_code), !anyNA(app_data$value))

saveRDS(app_data, file.path(OUT_DIR, "app_data.rds"))

################################################################################
# §4. Sanity checks
################################################################################

cat("=== Sanity checks ===\n")
n_expect <- 198 * (nrow(OUTCOMES) * 2 + nrow(COVARIATES) * 2 +
                   length(MODEL_IDS) * length(CRITERIA) * nrow(OUTCOMES) * 2)
cat(sprintf("  app_data rows: %d (expect %d)\n", nrow(app_data), n_expect))
stopifnot(nrow(app_data) == n_expect)

chk <- app_data %>%
  group_by(category) %>%
  summarise(n = n(), years = paste(sort(unique(year)), collapse = "/"),
            n_var = n_distinct(variable),
            n_model = n_distinct(model[!is.na(model)]),
            criteria = paste(sort(unique(criterion)), collapse = "/"),
            n_muni = n_distinct(muni_code), .groups = "drop")
print(as.data.frame(chk))

# Year asymmetry: residual/fitted are YEAR_RESID only
yr_resid <- app_data %>% filter(category %in% c("residual", "fitted")) %>%
  pull(year) %>% unique()
stopifnot(identical(yr_resid, YEAR_RESID))
cat(sprintf("  OK: residual/fitted are %d-only\n", YEAR_RESID))

# One label per (category, variable) -> selector choices do not duplicate
n_lab <- app_data %>% distinct(category, variable, label) %>%
  count(category, variable) %>% pull(n)
stopifnot(all(n_lab == 1))
cat("  OK: one label per (category, variable)\n")

# pAIC vs pBIC residuals: identical where gamma is the same in both criteria
# (1111, 0111, conventional models), different where it is not.
cat("  pAIC vs pBIC residuals, max |diff| over both outcomes:\n")
EXPECT_DIFF <- c("1011", "1101", "0101", "1001")
res_a <- app_data %>% filter(category == "residual", criterion == "pAIC")
res_b <- app_data %>% filter(category == "residual", criterion == "pBIC")
res_ab <- inner_join(res_a, res_b, by = c("muni_code", "outcome", "model"),
                     suffix = c("_a", "_b"))
stopifnot(nrow(res_ab) == nrow(res_a))
dmax <- res_ab %>% group_by(model) %>%
  summarise(max_abs_diff = max(abs(value_a - value_b)), .groups = "drop")
dmax <- dmax[match(MODEL_IDS, dmax$model), ]
for (i in seq_len(nrow(dmax))) {
  cat(sprintf("    %-5s %-8s %.6f\n", dmax$model[i],
              if (dmax$model[i] %in% EXPECT_DIFF) "(differ)" else "(same)",
              dmax$max_abs_diff[i]))
}
stopifnot(all(dmax$max_abs_diff[dmax$model %in% EXPECT_DIFF] > 0),
          all(dmax$max_abs_diff[!dmax$model %in% EXPECT_DIFF] == 0))
cat("  OK: pAIC/pBIC residuals identical for 1111/0111/conventional, differ for 1011/1101/0101/1001\n")

# Response / covariate ranges (pooled over both years) vs definition doc
REF <- list(aging65 = c(-1.8887, 3.5265), secondary_ind = c(-2.1456, 2.9951),
            ln_taxable_income = c(7.7212, 8.7828), ln_foreign_pop = c(0, 8.6415),
            net_inmigration = c(-5.79, 1.41), commuter_inflow = c(7.8, 231.1))
for (v in names(REF)) {
  r <- range(app_data$value[app_data$variable == v &
                              app_data$category %in% c("response", "covariate")])
  stopifnot(all(abs(round(r, 4) - REF[[v]]) < 5e-5))
  cat(sprintf("  OK: %-18s range [%.4f, %.4f] matches definition doc\n", v, r[1], r[2]))
}

# Cross-check against the reused muni_sf.rds (codes must align)
msf <- readRDS(file.path(OUT_DIR, "muni_sf.rds"))
tgt_codes <- msf$muni_code[msf$is_target]
stopifnot(setequal(unique(app_data$muni_code), tgt_codes))
cat("  OK: app_data muni_codes == reused muni_sf target layer (198)\n")

cat("\nDone. (muni_sf.rds / pref_borders.rds were NOT modified.)\n")
