################################################################################
# constants.R
# Single source of truth for variables, years, estimation criteria and models.
# Sourced by BOTH the prep scripts (prep/*.R) and the app (app.R / helpers.R),
# so labels and keys cannot drift apart. Base R only (no package needed).
#
# Data: Data_I (Kansai 198 municipalities), definitions in
#   data/real/data_source_and_processing_en.md  (sec. 3.1)
# Estimation: penalized model selection (pAIC / pBIC), results in output/real/
################################################################################

# --- Years --------------------------------------------------------------------
# time index in the CSVs -> fiscal year shown in the UI
YEARS      <- c(`1` = 2015L, `2` = 2020L)
YEAR_BASE  <- YEARS[["1"]]
YEAR_RESID <- YEARS[["2"]]   # residuals / fitted values exist for t = 2 only

# --- Estimation criteria --------------------------------------------------------
# Each penalized model has its own optimal gamma per criterion, so fitted values
# and residuals can differ between pAIC and pBIC.
CRITERIA        <- c("pAIC", "pBIC")
CRITERION_LABEL <- c(pAIC = "pAIC (penalized AIC)", pBIC = "pBIC (penalized BIC)")

# --- Outcomes (response variables) ---------------------------------------------
#   outcome : column prefix in the CSVs (y1 / y2)
#   key     : app-internal variable key
#   label   : full label (legends, popups)
#   short   : compact label (radio buttons, plot titles)
#   palette : viridisLite option (matches reference/create_heatmaps_en.r)
OUTCOMES <- data.frame(
  outcome = c("y1", "y2"),
  key     = c("aging65", "secondary_ind"),
  label   = c("Population Aged 65+ (%, standardized)",
              "Secondary-Industry Workers (%, standardized)"),
  short   = c("Aged 65+", "Secondary industry"),
  palette = c("plasma", "mako"),
  stringsAsFactors = FALSE
)

# --- Covariates -------------------------------------------------------------------
COVARIATES <- data.frame(
  col   = c("x_common1", "x_common2", "x_specific1_1", "x_specific2_1"),
  key   = c("ln_taxable_income", "ln_foreign_pop", "net_inmigration",
            "commuter_inflow"),
  label = c("Log Taxable Income per Taxpayer",
            "Log(1 + Foreign Population per 100k)",
            "Net In-migration Rate, Japanese migrants (%)",
            "Commuters from Other Municipalities (%)"),
  stringsAsFactors = FALSE
)

# Lookups by outcome id (y1 / y2)
OUTCOME_KEY   <- setNames(OUTCOMES$key,   OUTCOMES$outcome)
OUTCOME_LABEL <- setNames(OUTCOMES$label, OUTCOMES$outcome)
OUTCOME_SHORT <- setNames(OUTCOMES$short, OUTCOMES$outcome)

# --- Variable categories (Explore tab) -------------------------------------------
CATEGORY_LABELS <- c(
  response      = "Response variable",
  covariate     = "Covariate",
  fitted        = "Fitted value",
  residual      = "Residual",
  residual_diff = "Residual difference"
)
# Categories that exist for YEAR_RESID only and depend on the criterion.
RESID_CATS <- c("residual", "fitted", "residual_diff")

# --- Models ------------------------------------------------------------------------
# All 11 estimated models. Type label for display.
MODEL_TYPE <- c(
  "1111" = "MGNS",  "0111" = "MSEM",  "1011" = "MSAR",
  "1101" = "MGNS",  "0101" = "MSEM",  "1001" = "MSAR",
  "d0dd" = "IndSAR", "0ddd" = "IndSEM", "dddd" = "IndGNS",
  "0011" = "VARX",  "000d" = "IndReg"
)
MODEL_IDS <- names(MODEL_TYPE)

MODEL_LABELS <- setNames(sprintf("%s (%s)", MODEL_TYPE, MODEL_IDS), MODEL_IDS)

# Comparison axis: new methodology vs conventional individual regression.
MODEL_FAMILY <- c(
  "1111" = "new", "0111" = "new", "1011" = "new",
  "1101" = "new", "0101" = "new", "1001" = "new",
  "d0dd" = "conventional", "0ddd" = "conventional", "dddd" = "conventional",
  "0011" = "conventional", "000d" = "conventional"
)
FAMILY_LABEL <- c(new = "New methodology", conventional = "Conventional (individual)")
FAMILY_COLOR <- c(new = "#1b7837", conventional = "#762a83")  # green=new, purple=conv.
