################################################################################
# build_spatial_diag.R
#
# Pre-compute the spatial-dependence diagnostics consumed by the Shiny app:
#   - B-1: global Moran's I of each model's residuals (permutation p-value)
#   - B-3: LISA (local Moran) quadrant + permutation p-value per municipality
# for every model x outcome x criterion (pAIC / pBIC).
#
# Residuals are the INNOVATION residuals `y*_resid` (never `y*_resid_reduced`).
#
# Design principle: the APP READS ONLY. All spatial statistics are computed
# here and stored in data/spatial_diag.rds. The app never runs spdep.
#
# Inputs:
#   data/real/W_row_standardized_198_for_R.csv   row-standardized W (198x198,
#                                                 row order = region 1..198)
#   data/real/region_id_mapping198.csv           region, code (6-digit), name
#   output/real/fitted_residuals_<criterion>/fitted_residuals_<ID>.csv
#
# Output:
#   shiny_app/data/spatial_diag.rds  = list(moran = <df>, lisa = <df>)
#
# Run from the repository root:
#   Rscript shiny_app/prep/build_spatial_diag.R
################################################################################

suppressMessages({
  library(dplyr)
  library(spdep)
})

set.seed(20260614)   # reproducible permutation p-values
NSIM <- 999

REPO_ROOT <- Sys.getenv("MSNR_REPO_ROOT", getwd())
source(file.path(REPO_ROOT, "shiny_app", "R", "constants.R"))

PATH_W       <- file.path(REPO_ROOT, "data", "real", "W_row_standardized_198_for_R.csv")
PATH_MAPPING <- file.path(REPO_ROOT, "data", "real", "region_id_mapping198.csv")
RESID_DIR    <- function(crit) {
  file.path(REPO_ROOT, "output", "real", paste0("fitted_residuals_", crit))
}
OUT_DIR      <- file.path(REPO_ROOT, "shiny_app", "data")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

cat("=== build_spatial_diag.R (Moran + LISA; app reads only) ===\n")
cat("  REPO_ROOT:", REPO_ROOT, "\n")
stopifnot(file.exists(PATH_W), file.exists(PATH_MAPPING),
          all(dir.exists(vapply(CRITERIA, RESID_DIR, character(1)))))

################################################################################
# §1. region -> code_5 and the spatial weights listw
################################################################################

mapping <- read.csv(PATH_MAPPING, fileEncoding = "UTF-8")
stopifnot(identical(names(mapping), c("region", "code", "name")))
mapping$code   <- as.character(mapping$code)
mapping$code_5 <- substr(mapping$code, 1, 5)
stopifnot(nrow(mapping) == 198)
reg2code <- setNames(mapping$code_5, mapping$region)

# W rows are in region order 1..198. Build a listw; the matrix is already
# row-standardized, style="W" is idempotent here.
W <- as.matrix(read.csv(PATH_W))
stopifnot(nrow(W) == 198, ncol(W) == 198)
listw <- mat2listw(W, style = "W", zero.policy = TRUE)

# code_5 aligned to region order 1..198 (= W row order).
code_by_region <- unname(reg2code[as.character(1:198)])
stopifnot(!any(is.na(code_by_region)))

################################################################################
# §2. Helpers: read a model's residual vector aligned to region 1..198
################################################################################

resid_aligned <- function(crit, mid, oc) {
  fpath <- file.path(RESID_DIR(crit), sprintf("fitted_residuals_%s.csv", mid))
  if (!file.exists(fpath)) stop(sprintf("Missing: %s", fpath))
  df <- read.csv(fpath)
  df <- df[df$time == 2, , drop = FALSE]                 # YEAR_RESID only
  stopifnot(nrow(df) == 198)
  col <- paste0(oc, "_resid")                            # innovation residual
  x <- df[[col]][match(1:198, df$region)]                # aligned to W row order
  stopifnot(!anyNA(x))
  x
}

# Pick the simulated p-value column from a localmoran_perm result robustly.
lisa_pcol <- function(lmp) {
  nm <- colnames(lmp)
  if ("Pr(folded) Sim" %in% nm)       return("Pr(folded) Sim")
  if ("Pr(z != E(Ii)) Sim" %in% nm)   return("Pr(z != E(Ii)) Sim")
  sim <- grep("Sim", nm, value = TRUE)
  if (length(sim)) return(tail(sim, 1))
  tail(nm, 1)
}

################################################################################
# §3. Global Moran's I (B-1) and LISA (B-3) for each criterion x model x outcome
################################################################################

moran_rows <- list()
lisa_rows  <- list()

for (crit in CRITERIA) {
  for (mid in MODEL_IDS) {
    for (oc in OUTCOMES$outcome) {
      x <- resid_aligned(crit, mid, oc)
      id <- paste(crit, mid, oc)

      # --- B-1: global Moran's I (analytical I/E/Var/z) + permutation p ---
      mt <- moran.test(x, listw, zero.policy = TRUE, randomisation = TRUE)
      mc <- moran.mc(x, listw, nsim = NSIM, zero.policy = TRUE)
      moran_rows[[id]] <- data.frame(
        model       = mid,
        model_type  = unname(MODEL_TYPE[mid]),
        family      = unname(MODEL_FAMILY[mid]),
        outcome     = oc,
        criterion   = crit,
        morans_i    = unname(mt$estimate["Moran I statistic"]),
        expectation = unname(mt$estimate["Expectation"]),
        variance    = unname(mt$estimate["Variance"]),
        z           = unname(mt$statistic),
        p_analytic  = unname(mt$p.value),
        p_value     = unname(mc$p.value),     # permutation (preferred)
        n_sim       = NSIM,
        stringsAsFactors = FALSE
      )

      # --- B-3: LISA quadrant + permutation p-value (alpha applied in app) ---
      lmp  <- localmoran_perm(x, listw, nsim = NSIM, zero.policy = TRUE)
      pcol <- lisa_pcol(lmp)
      pval <- as.numeric(lmp[, pcol])
      lagx <- lag.listw(listw, x, zero.policy = TRUE)
      xc   <- x - mean(x, na.rm = TRUE)
      lc   <- lagx - mean(lagx, na.rm = TRUE)
      quad <- ifelse(xc > 0 & lc > 0, "HH",
               ifelse(xc < 0 & lc < 0, "LL",
                ifelse(xc > 0 & lc < 0, "HL", "LH")))

      lisa_rows[[id]] <- data.frame(
        muni_code = code_by_region,
        model     = mid,
        outcome   = oc,
        criterion = crit,
        local_i   = as.numeric(lmp[, "Ii"]),
        quadrant  = quad,      # HH/LL/HL/LH (raw); NS derived in app via p<alpha
        p_value   = pval,
        stringsAsFactors = FALSE
      )
    }
  }
}

moran <- bind_rows(moran_rows)
lisa  <- bind_rows(lisa_rows)

saveRDS(list(moran = moran, lisa = lisa),
        file.path(OUT_DIR, "spatial_diag.rds"))

################################################################################
# §4. Sanity checks
################################################################################

cat("=== Sanity checks ===\n")
n_m <- length(MODEL_IDS) * nrow(OUTCOMES) * length(CRITERIA)
cat(sprintf("  moran rows: %d (expect %d = 11 models x 2 outcomes x 2 criteria)\n",
            nrow(moran), n_m))
cat(sprintf("  lisa  rows: %d (expect %d = %d x 198)\n", nrow(lisa), n_m * 198, n_m))
stopifnot(nrow(moran) == n_m, nrow(lisa) == n_m * 198)
stopifnot(all(!is.na(lisa$muni_code)), all(is.finite(lisa$p_value)))

# Reference: reference/heatmaps_Data_I_v2_console.txt (pAIC, moran.test I)
REF_I <- data.frame(model = c("1111", "1111", "dddd", "dddd"),
                    outcome = c("y1", "y2", "y1", "y2"),
                    ref = c(0.0169, 0.0223, 0.0073, 0.0110))
cat("  Moran's I vs reference console (pAIC, 4 decimals):\n")
for (i in seq_len(nrow(REF_I))) {
  got <- moran$morans_i[moran$criterion == "pAIC" & moran$model == REF_I$model[i] &
                          moran$outcome == REF_I$outcome[i]]
  ok <- round(got, 4) == REF_I$ref[i]
  cat(sprintf("    e%s(%s)  I = %+.6f  -> %+.4f  ref %+.4f  %s\n",
              substr(REF_I$outcome[i], 2, 2), REF_I$model[i], got, round(got, 4),
              REF_I$ref[i], if (ok) "OK" else "NG"))
  stopifnot(ok)
}

cat("  Moran's I by family x criterion (mean |I|, n with perm p < 0.05):\n")
print(moran %>% group_by(criterion, family, outcome) %>%
        summarise(mean_absI = mean(abs(morans_i)),
                  n_sig = sum(p_value < 0.05), n = n(), .groups = "drop") %>%
        as.data.frame())

cat("\nDone. -> shiny_app/data/spatial_diag.rds\n")
