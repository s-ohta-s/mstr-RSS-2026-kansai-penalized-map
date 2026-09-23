################################################################################
# helpers.R
# Label generation, color-scale selection, hover labels, popup tables and the
# ggplot builders. Sourced by app.R AFTER constants.R (variables, years,
# criteria and model tables live there). Pure functions; no global state
# except the palette cache.
#
# Penalized update: residual/fitted rows carry a `criterion` (pAIC / pBIC).
# Residual differences are computed on the fly for any model pair within one
# criterion.
################################################################################

suppressMessages({
  library(dplyr)
  library(htmltools)
  library(leaflet)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

# Pre-generated viridis color ramps (passed to leaflet::colorNumeric as vectors
# so palette names need not be recognized by the installed leaflet version).
.PAL <- list(
  plasma  = viridisLite::plasma(256),
  mako    = viridisLite::mako(256),
  viridis = viridisLite::viridis(256)
)
NA_COLOR <- "#d9d9d9"

# Outcome id (y1 / y2) -> color ramp, from constants.R OUTCOMES$palette.
outcome_pal <- function(outcome) .PAL[[OUTCOMES$palette[OUTCOMES$outcome == outcome]]]

# Named choices for outcome radio buttons (labels = short names, values = y1/y2).
OUTCOME_CHOICES   <- setNames(OUTCOMES$outcome, OUTCOMES$short)
CRITERION_CHOICES <- setNames(CRITERIA, CRITERION_LABEL[CRITERIA])

#' Format a numeric value for display.
fmt_value <- function(x) {
  ifelse(is.na(x), "–", formatC(x, format = "f", digits = 3))
}

#' Variable choices (named vector: names = labels shown, values = variable keys)
#' for a given category, in data order. Residual/fitted labels carry no
#' criterion, so pAIC and pBIC rows collapse to one choice per variable.
var_choices <- function(app_data, category) {
  v <- app_data %>%
    filter(.data$category == !!category) %>%
    distinct(.data$variable, .data$label)
  setNames(v$variable, v$label)
}

#' Model choices (named vector: names = labels, values = model ids) restricted to
#' the models actually present in the residual rows of app_data, in MODEL_IDS
#' order.
model_choices <- function(app_data) {
  present <- unique(app_data$model[app_data$category == "residual"])
  ids <- MODEL_IDS[MODEL_IDS %in% present]
  setNames(ids, unname(MODEL_LABELS[ids]))
}

#' Pre-compute shared color-scale domains from the full long data.
#'  - response & fitted share one range per outcome (both criteria included);
#'  - covariates use a per-variable viridis range;
#'  - residuals use one symmetric RdBu range across all outcomes/models/criteria.
#' Residual-difference scales are computed dynamically (per selected pair) in
#' the server, so they are not part of this object.
compute_scales <- function(app_data) {
  rf <- app_data %>% filter(.data$category %in% c("response", "fitted"))
  resp <- lapply(setNames(OUTCOMES$outcome, OUTCOMES$outcome), function(oc)
    range(rf$value[rf$outcome == oc], na.rm = TRUE))

  covs <- app_data %>%
    filter(.data$category == "covariate") %>%
    group_by(.data$variable) %>%
    summarise(lo = min(value, na.rm = TRUE), hi = max(value, na.rm = TRUE),
              .groups = "drop")
  cov_ranges <- setNames(
    lapply(seq_len(nrow(covs)), function(i) c(covs$lo[i], covs$hi[i])),
    covs$variable
  )

  rmax <- max(abs(app_data$value[app_data$category == "residual"]), na.rm = TRUE)

  list(resp = resp, cov = cov_ranges, resid = c(-rmax, rmax))
}

#' Pick domain + palette for the current selection.
#' (residual_diff is handled directly in the server with a dynamic domain.)
scale_for <- function(scales, category, variable, outcome) {
  if (category == "residual") {
    list(domain = scales$resid, palette = "RdBu", diverging = TRUE)
  } else if (category == "covariate") {
    dom <- scales$cov[[variable]] %||% c(0, 1)
    list(domain = dom, palette = .PAL$viridis, diverging = FALSE)
  } else { # response / fitted
    list(domain = scales$resp[[outcome]], palette = outcome_pal(outcome),
         diverging = FALSE)
  }
}

#' Symmetric diverging scale spec for a residual-difference vector.
diff_scale <- function(values) {
  m <- max(abs(values), na.rm = TRUE)
  if (!is.finite(m) || m == 0) m <- 1
  list(domain = c(-m, m), palette = "RdBu", diverging = TRUE)
}

#' Build a leaflet colorNumeric function for a scale spec.
#' For diverging scales the brewer "RdBu" ramp maps low -> red, high -> blue,
#' matching create_heatmaps_en.r (scale_fill_distiller(palette="RdBu", direction=1)).
make_pal <- function(sp) {
  colorNumeric(palette = sp$palette, domain = sp$domain, na.color = NA_COLOR)
}

#' Legend-only palette with the color ramp reversed. Combined with
#' labelFormat(transform = sort decreasing) in addLegend, this flips BOTH the
#' gradient and the tick labels together, so small (dark) values sit at the
#' bottom while color<->value association stays correct. The map polygons keep
#' using make_pal().
make_legend_pal <- function(sp) {
  cols <- if (is.character(sp$palette) && length(sp$palette) == 1) {
    rev(RColorBrewer::brewer.pal(11, sp$palette))  # brewer name, e.g. "RdBu"
  } else {
    rev(sp$palette)                                # explicit color vector
  }
  colorNumeric(palette = cols, domain = sp$domain, na.color = NA_COLOR)
}

#' Hover labels for the target layer (one HTML string per municipality, aligned
#' to `target` row order). The list must be UNNAMED: a named list is serialized
#' to a JSON object, and leaflet then hands every tooltip the whole object.
hover_labels <- function(target, values, var_label) {
  mapply(function(nm, ja, v) {
    HTML(sprintf("<b>%s</b> <span style='color:#888'>(%s)</span><br/>%s: %s",
                 nm, ja, var_label, fmt_value(v)))
  }, target$muni_name, target$muni_name_ja, values,
  SIMPLIFY = FALSE, USE.NAMES = FALSE)
}

#' Detailed click popups for ALL target municipalities at once (vector aligned
#' to `target` row order): response/covariates at the selected year, plus the
#' YEAR_RESID per-model residuals and fitted values for the selected criterion.
#' Built column-wise (one lookup per variable), so it stays fast on redraw.
popup_all <- function(app_data, target, year, criterion) {
  codes <- target$muni_code
  pick <- function(cat, var, yr, crit = NA_character_) {
    keep <- app_data$category == cat & app_data$variable == var &
      app_data$year == yr
    if (!is.na(crit)) keep <- keep & app_data$criterion %in% crit
    r <- app_data[keep, ]
    r$value[match(codes, r$muni_code)]
  }
  trow <- function(lab, val) {
    sprintf("<tr><td style='padding-right:10px'>%s</td><td style='text-align:right'>%s</td></tr>",
            lab, fmt_value(val))
  }
  trow3 <- function(lab, v1, v2) {
    sprintf(paste0("<tr><td style='padding-left:14px;padding-right:10px'>%s</td>",
                   "<td style='text-align:right;padding-right:8px'>%s</td>",
                   "<td style='text-align:right'>%s</td></tr>"),
            lab, fmt_value(v1), fmt_value(v2))
  }
  section <- function(title, span = 2) {
    sprintf("<tr><td colspan='%d' style='padding-top:6px;font-weight:bold;color:#444'>%s</td></tr>",
            span, title)
  }
  oc <- OUTCOMES$outcome  # exactly two outcomes -> two value columns
  per_model <- function(cat) {
    lapply(MODEL_IDS, function(mid) {
      trow3(MODEL_LABELS[[mid]],
            pick(cat, paste0(OUTCOME_KEY[[oc[1]]], "_", mid), YEAR_RESID, criterion),
            pick(cat, paste0(OUTCOME_KEY[[oc[2]]], "_", mid), YEAR_RESID, criterion))
    })
  }
  oc_head <- paste(OUTCOME_SHORT[oc], collapse = " | ")

  pieces <- c(
    list(section(sprintf("Response (%d)", year))),
    lapply(seq_len(nrow(OUTCOMES)), function(i)
      trow(OUTCOMES$label[i], pick("response", OUTCOMES$key[i], year))),
    list(section(sprintf("Covariates (%d)", year))),
    lapply(seq_len(nrow(COVARIATES)), function(i)
      trow(COVARIATES$label[i], pick("covariate", COVARIATES$key[i], year))),
    list(section(sprintf("Residuals %d, %s — %s", YEAR_RESID, criterion, oc_head), span = 3)),
    per_model("residual"),
    list(section(sprintf("Fitted %d, %s — %s", YEAR_RESID, criterion, oc_head), span = 3)),
    per_model("fitted")
  )
  rows <- do.call(paste0, pieces)

  sprintf(
    "<div style='font-size:12px;line-height:1.35'><div style='font-size:14px;font-weight:bold;margin-bottom:2px'>%s <span style='color:#888;font-weight:normal'>(%s)</span></div><table>%s</table></div>",
    target$muni_name, target$muni_name_ja, rows
  )
}

################################################################################
# Model family, LISA cluster derivation, and ggplot builders for B-1 (Moran)
# and C-1 (estimated matrices). All pure functions; no estimation.
################################################################################

# --- Estimated-matrix display names (C-1) -------------------------------------
MATRIX_LABELS <- c(
  A      = "A — temporal autoregression",
  R      = "R — spatial lag (SAR) coef.",
  Lambda = "Λ — spatial error (SEM) coef.",
  Sigma  = "Σ — residual covariance"
)

# --- LISA cluster colors (GeoDa convention) ----------------------------------
LISA_LEVELS <- c("HH", "LL", "HL", "LH", "NS")
LISA_COLORS <- c(HH = "#b2182b", LL = "#2166ac", HL = "#ef8a62",
                 LH = "#67a9cf", NS = "#e6e6e6")
LISA_LABELS <- c(HH = "High-High", LL = "Low-Low", HL = "High-Low",
                 LH = "Low-High", NS = "Not significant")

#' Model choices as an optgroup list (Conventional first, then New methodology),
#' restricted to models present in app_data residual rows.
model_choices_grouped <- function(app_data) {
  ids <- unname(model_choices(app_data))
  fam <- MODEL_FAMILY[ids]
  mk <- function(f) {
    sub <- ids[fam == f]
    setNames(as.list(sub), unname(MODEL_LABELS[sub]))
  }
  out <- list()
  if (any(fam == "conventional")) out[[FAMILY_LABEL["conventional"]]] <- mk("conventional")
  if (any(fam == "new"))          out[[FAMILY_LABEL["new"]]]          <- mk("new")
  out
}

#' Derive the final LISA cluster factor from the stored quadrant + p-value at a
#' chosen significance level (NS when not significant). Returns an ordered factor.
lisa_cluster <- function(quadrant, p_value, alpha = 0.05) {
  cl <- ifelse(is.finite(p_value) & p_value < alpha, quadrant, "NS")
  factor(cl, levels = LISA_LEVELS)
}

#' Symmetric diverging palette for the Panel-change (Delta) map.
delta_scale <- function(values) {
  m <- max(abs(values), na.rm = TRUE)
  if (!is.finite(m) || m == 0) m <- 1
  list(domain = c(-m, m), palette = "RdBu", diverging = TRUE)
}

# ---------------------------------------------------------------------------
# B-1: global Moran's I bar chart (residual spatial autocorrelation by model)
# ---------------------------------------------------------------------------
build_moran_plot <- function(moran, outcome, criterion) {
  d <- moran %>% filter(.data$outcome == !!outcome, .data$criterion == !!criterion)
  d$label <- MODEL_LABELS[d$model]
  d$sig   <- ifelse(d$p_value < 0.05, "*", "")
  d <- d[order(d$morans_i), ]
  d$label <- factor(d$label, levels = d$label)
  ggplot(d, aes(x = .data$morans_i, y = .data$label, fill = .data$family)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, color = "#333333", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.3f%s", .data$morans_i, .data$sig),
                  hjust = ifelse(.data$morans_i >= 0, -0.12, 1.12)),
              size = 3.3, color = "#222222") +
    scale_fill_manual(values = FAMILY_COLOR, labels = FAMILY_LABEL,
                      name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
    labs(x = "Moran's I of residuals  (* p < 0.05, permutation)", y = NULL,
         title = sprintf("Residual spatial autocorrelation — %s (%d)",
                         OUTCOME_SHORT[[outcome]], YEAR_RESID),
         subtitle = sprintf("%s. Closer to 0 = more spatial structure removed",
                            CRITERION_LABEL[[criterion]])) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "top",
          panel.grid.major.y = element_blank(),
          plot.title = element_text(face = "bold"))
}

# ---------------------------------------------------------------------------
# C-1: estimated-matrix heatmaps (A / R / Lambda / Sigma) for selected models
# ---------------------------------------------------------------------------
build_matrix_plot <- function(coef, models, criterion) {
  d <- coef %>% filter(.data$model %in% models, .data$criterion == !!criterion)
  d$model_lab  <- factor(MODEL_LABELS[d$model], levels = MODEL_LABELS[models])
  d$matrix_lab <- factor(MATRIX_LABELS[d$matrix], levels = MATRIX_LABELS)
  d$fill_val   <- ifelse(d$present, d$estimate, NA_real_)
  d$cell_text  <- ifelse(d$present,
                         paste0(formatC(d$estimate, format = "f", digits = 2),
                                d$signif),
                         "0")
  d$text_col   <- ifelse(!is.na(d$fill_val) & abs(d$fill_val) >
                           0.6 * max(abs(d$fill_val), na.rm = TRUE),
                         "white", "#333333")
  mabs <- max(abs(d$fill_val), na.rm = TRUE); if (!is.finite(mabs)) mabs <- 1

  ggplot(d, aes(x = factor(.data$col), y = factor(.data$row))) +
    geom_tile(aes(fill = .data$fill_val), color = "white", linewidth = 1) +
    geom_text(aes(label = .data$cell_text, color = I(.data$text_col)), size = 3.6) +
    scale_y_discrete(limits = c("2", "1")) +              # row 1 on top
    scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac",
                         midpoint = 0, limits = c(-mabs, mabs),
                         na.value = "#dcdcdc", name = "estimate") +
    facet_grid(model_lab ~ matrix_lab, switch = "y") +
    labs(x = NULL, y = NULL,
         title = sprintf("Estimated parameter matrices — %s",
                         CRITERION_LABEL[[criterion]]),
         subtitle = paste0(
           "Grey = structurally 0 (diagonal constraint). */**/*** significance.\n",
           "Significance: Ψ-based SE; Hessian-based SE where Ψ is unavailable (R, Λ).\n",
           "No mark where no SE is reported (Σ, conventional R/Λ).")) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          strip.text.y.left = element_text(angle = 0, face = "bold"),
          strip.text.x = element_text(face = "bold", size = 9),
          axis.text = element_text(color = "#888888"),
          plot.title = element_text(face = "bold"),
          legend.position = "right")
}
