################################################################################
# app.R  —  MGNST result viewer (Kansai 198 municipalities), penalized (pAIC/pBIC)
#
# Theme: "Visual Comparison of New Spatial Methodology with Conventional
#         Individual Regression Models."
#
# Reads pre-computed artifacts only (built by prep/*.R); performs NO estimation
# and NO spatial statistics at runtime.
#
#   R/constants.R           variables, years, criteria, model tables
#   data/muni_sf.rds        two-layer polygons (198 target + background prefs)
#   data/pref_borders.rds   prefecture outlines
#   data/app_data.rds       long response/covariate/residual/fitted (+ criterion)
#   data/spatial_diag.rds   list(moran, lisa)         — B-1 / B-3 (built in prep)
#   data/coef_matrices.rds  list(coef, fit, gamma)    — C-1       (built in prep)
#
# Residuals are the innovation residuals of the penalized estimation.
# Fitted values and residuals exist for each estimation criterion (pAIC / pBIC);
# every tab has a criterion selector, and the swipe tab picks one per side.
#
# Tabs:
#   1 Explore              single map; response / covariate / fitted / residual /
#                          residual difference
#   2 Panel change         Delta = 2020 - 2015 map (response / covariate)
#   3 Compare (swipe)      side-by-side residual/fitted of two model x criterion
#   4 Spatial diagnostics  Moran's I bar (B-1) + LISA cluster map (B-3)
#   5 Model structure      estimated A/R/Lambda/Sigma matrices (C-1) + fit table
#
# Run from the repository root:  shiny::runApp("shiny_app")
################################################################################

suppressMessages({
  library(shiny)
  library(leaflet)
  library(leaflet.extras2)
  library(sf)
  library(dplyr)
  library(htmltools)
  library(ggplot2)
})

# --- Load data and helpers ONCE (shared across all sessions) -----------------
app_dir <- if (file.exists("app.R")) "." else "shiny_app"
source(file.path(app_dir, "R", "constants.R"))
source(file.path(app_dir, "R", "helpers.R"))

muni_sf      <- readRDS(file.path(app_dir, "data", "muni_sf.rds"))
pref_borders <- readRDS(file.path(app_dir, "data", "pref_borders.rds"))
app_data     <- readRDS(file.path(app_dir, "data", "app_data.rds"))
spatial_diag <- readRDS(file.path(app_dir, "data", "spatial_diag.rds"))
coef_mat     <- readRDS(file.path(app_dir, "data", "coef_matrices.rds"))

moran_df <- spatial_diag$moran
lisa_df  <- spatial_diag$lisa
coef_df  <- coef_mat$coef
fit_df   <- coef_mat$fit

target     <- muni_sf[muni_sf$is_target, ]
background <- muni_sf[!muni_sf$is_target, ]
scales     <- compute_scales(app_data)
MODEL_CH     <- model_choices(app_data)          # flat (residual_diff selectors)
MODEL_CH_GRP <- model_choices_grouped(app_data)  # optgroup (swipe/diag/structure)

KANSAI_PREF_CODES <- c("25", "26", "27", "28", "29", "30")
kansai_borders <- pref_borders[pref_borders$pref_code %in% KANSAI_PREF_CODES, ]
bb <- as.numeric(st_bbox(target))

MAP_H <- "calc(100vh - 130px)"  # leaflet height inside a navbar tab

# --- Basemap -------------------------------------------------------------------
# CARTO Positron. With CARTO_BASEMAP_KEY set (local .Renviron, or a secret on
# Posit Connect Cloud) the keyed CARTO URL is used; without it the app falls
# back to leaflet's built-in CartoDB.Positron provider, so a fresh clone still
# runs. The key itself is never printed.
CARTO_BASEMAP_KEY <- trimws(Sys.getenv("CARTO_BASEMAP_KEY"))
HAS_CARTO_KEY     <- nzchar(CARTO_BASEMAP_KEY)

CARTO_POSITRON_URL <- paste0(
  "https://{s}.basemaps.cartocdn.com/light_all/",
  "{z}/{x}/{y}{r}.png?key=",
  utils::URLencode(CARTO_BASEMAP_KEY, reserved = TRUE)
)
CARTO_ATTRIBUTION <- paste0(
  '&copy; <a href="https://www.openstreetmap.org/copyright">',
  'OpenStreetMap</a> contributors, ',
  '&copy; <a href="https://carto.com/attributions">CARTO</a>'
)

message(if (HAS_CARTO_KEY) "[basemap] CARTO Positron (API key from CARTO_BASEMAP_KEY)"
        else "[basemap] CARTO_BASEMAP_KEY not set -> fallback providers$CartoDB.Positron")

#' Add the basemap tile layer. `layerId` / `pane` are used by the swipe tab:
#' leaflet.extras2 side-by-side looks its layers up as TILE layers by layerId
#' and clips by pane, so both code paths must register a tile layer.
add_basemap <- function(map, layerId = NULL, pane = NULL) {
  if (HAS_CARTO_KEY) {
    opts <- tileOptions(subdomains = "abcd", maxZoom = 20)
    if (!is.null(pane)) opts$pane <- pane
    addTiles(map, urlTemplate = CARTO_POSITRON_URL,
             attribution = CARTO_ATTRIBUTION, layerId = layerId, options = opts)
  } else {
    opts <- providerTileOptions(maxZoom = 20)
    if (!is.null(pane)) opts$pane <- pane
    addProviderTiles(map, providers$CartoDB.Positron, layerId = layerId,
                     options = opts)
  }
}

#' Base leaflet map shared by the single-map tabs.
base_map <- function() {
  leaflet(options = leafletOptions(minZoom = 6)) %>%
    add_basemap() %>%
    addPolylines(data = kansai_borders, color = "#888888",
                 weight = 0.8, opacity = 0.6) %>%
    fitBounds(bb[1], bb[2], bb[3], bb[4])
}

# --- Small data accessors (aligned to `target` row order) --------------------

# residual/fitted value vector for one outcome key/model/criterion, YEAR_RESID.
vec_for <- function(categ, oc_key, mid, crit) {
  r <- app_data %>%
    filter(.data$category == categ,
           .data$variable == paste0(oc_key, "_", mid),
           .data$criterion == crit,
           .data$year == YEAR_RESID)
  r$value[match(target$muni_code, r$muni_code)]
}
resid_vec <- function(oc_key, mid, crit) vec_for("residual", oc_key, mid, crit)

# Delta (YEAR_RESID - YEAR_BASE) for a response/covariate variable, aligned.
delta_vec <- function(var) {
  d1 <- app_data %>% filter(.data$variable == var, .data$year == YEAR_RESID)
  d0 <- app_data %>% filter(.data$variable == var, .data$year == YEAR_BASE)
  v1 <- d1$value[match(target$muni_code, d1$muni_code)]
  v0 <- d0$value[match(target$muni_code, d0$muni_code)]
  list(d = v1 - v0, v0 = v0, v1 = v1)
}

# Global Moran's I lookup (for the swipe annotation).
get_moran <- function(mid, oc, crit) {
  r <- moran_df %>% filter(.data$model == mid, .data$outcome == oc,
                           .data$criterion == crit)
  if (nrow(r) == 0) NA_real_ else r$morans_i[1]
}

PC_CH <- c(var_choices(app_data, "response"), var_choices(app_data, "covariate"))

################################################################################
# UI
################################################################################

busy_overlay <- tagList(
  tags$head(tags$style(HTML("
    #busy-overlay { display:none; position:fixed; top:0; left:0; right:0;
      bottom:0; background:rgba(255,255,255,0.5); z-index:10000;
      align-items:center; justify-content:center; pointer-events:none; }
    #busy-overlay.active { display:flex; }
    #busy-overlay .box { background:#fff; border:1px solid #ccc;
      border-radius:8px; padding:14px 20px; font-size:14px; color:#333;
      box-shadow:0 2px 10px rgba(0,0,0,.18); display:flex; align-items:center;
      gap:10px; }
    #busy-overlay .spinner { width:18px; height:18px; border:3px solid #ddd;
      border-top-color:#3c8dbc; border-radius:50%;
      animation:busy-spin .8s linear infinite; }
    @keyframes busy-spin { to { transform:rotate(360deg); } }
    .fam-badge { padding:2px 8px; border-radius:10px; color:#fff;
      font-weight:bold; font-size:12px; }
    .crit-badge { padding:1px 6px; border-radius:8px; background:#eee;
      color:#333; font-size:11px; font-weight:bold; }
  ")),
  tags$script(HTML("
    $(document).on('shiny:busy', function(){ $('#busy-overlay').addClass('active'); });
    $(document).on('shiny:idle', function(){ $('#busy-overlay').removeClass('active'); });
  "))),
  div(id = "busy-overlay",
      div(class = "box", div(class = "spinner"), span("Loading / rendering…")))
)

RESID_CATS_JS <- paste(RESID_CATS, collapse = "','")

# ---- Tab 1: Explore ---------------------------------------------------------
tab_explore <- tabPanel(
  "Explore",
  sidebarLayout(
    sidebarPanel(
      width = 3,
      radioButtons("exp_category", "Variable category",
        choices  = setNames(names(CATEGORY_LABELS), unname(CATEGORY_LABELS)),
        selected = "response"),
      conditionalPanel(
        sprintf("['%s'].indexOf(input.exp_category) >= 0", RESID_CATS_JS),
        radioButtons("exp_criterion", "Estimation criterion",
          choices = CRITERION_CHOICES, selected = "pAIC")),
      conditionalPanel(
        "input.exp_category != 'residual_diff'",
        selectInput("exp_variable", "Variable",
                    choices = var_choices(app_data, "response"))),
      conditionalPanel(
        "input.exp_category == 'residual_diff'",
        radioButtons("exp_diff_outcome", "Outcome",
          choices = OUTCOME_CHOICES, selected = "y1", inline = TRUE),
        selectInput("exp_diff_a", "Model A (minuend)", choices = MODEL_CH,
          selected = if ("0ddd" %in% MODEL_CH) "0ddd" else MODEL_CH[1]),
        selectInput("exp_diff_b", "Model B (subtrahend)", choices = MODEL_CH,
          selected = if ("0111" %in% MODEL_CH) "0111" else MODEL_CH[2]),
        helpText("Mapped value = e_A − e_B (same criterion).")),
      conditionalPanel(
        sprintf("['%s'].indexOf(input.exp_category) < 0", RESID_CATS_JS),
        radioButtons("exp_year", "Year", choices = unname(YEARS),
                     selected = YEAR_RESID, inline = TRUE)),
      conditionalPanel(
        sprintf("['%s'].indexOf(input.exp_category) >= 0", RESID_CATS_JS),
        helpText(sprintf("Residuals / fitted / differences are %d-only.",
                         YEAR_RESID))),
      tags$hr(),
      checkboxInput("exp_show_choropleth", "Show choropleth layer", TRUE),
      checkboxInput("exp_show_legend", "Show legend", TRUE),
      sliderInput("exp_opacity", "Choropleth opacity", 0, 1, 0.85, 0.05),
      helpText(tags$b("Hover"), "= name + value;", tags$b("click"), "= full table.")
    ),
    mainPanel(width = 9, leafletOutput("map", height = MAP_H))
  )
)

# ---- Tab 2: Panel change (Delta) -------------------------------------------
tab_change <- tabPanel(
  "Panel change",
  sidebarLayout(
    sidebarPanel(
      width = 3,
      helpText(tags$b(sprintf("Δ = value(%d) − value(%d)", YEAR_RESID, YEAR_BASE)),
               "for response & covariate variables. Diverging RdBu, 0-centered."),
      selectInput("pc_variable", "Variable", choices = PC_CH,
                  selected = PC_CH[1]),
      tags$hr(),
      checkboxInput("pc_show_legend", "Show legend", TRUE),
      sliderInput("pc_opacity", "Choropleth opacity", 0, 1, 0.85, 0.05),
      helpText(sprintf("Hover shows %d, %d and Δ.", YEAR_BASE, YEAR_RESID))
    ),
    mainPanel(width = 9, leafletOutput("map_delta", height = MAP_H))
  )
)

# ---- Tab 3: Compare (swipe) -------------------------------------------------
tab_swipe <- tabPanel(
  "Compare (swipe)",
  sidebarLayout(
    sidebarPanel(
      width = 3,
      helpText(tags$b("Drag the central handle"), "to wipe between two maps.",
               "Both sides share one color scale."),
      radioButtons("sw_outcome", "Outcome",
        choices = OUTCOME_CHOICES, selected = "y1", inline = TRUE),
      radioButtons("sw_quantity", "Quantity",
        choices = c("Residual" = "residual", "Fitted" = "fitted"),
        selected = "residual", inline = TRUE),
      tags$hr(),
      selectInput("sw_left", "Left model", choices = MODEL_CH_GRP,
        selected = "0ddd"),
      radioButtons("sw_left_crit", "Left criterion", choices = CRITERIA,
        selected = "pAIC", inline = TRUE),
      tags$hr(),
      selectInput("sw_right", "Right model", choices = MODEL_CH_GRP,
        selected = "0111"),
      radioButtons("sw_right_crit", "Right criterion", choices = CRITERIA,
        selected = "pAIC", inline = TRUE),
      tags$hr(),
      checkboxInput("sw_show_legend", "Show legend", TRUE),
      sliderInput("sw_opacity", "Choropleth opacity", 0, 1, 0.9, 0.05),
      helpText("Default: conventional 0ddd (left) vs new 0111 (right), both pAIC.",
               "Pick the same model with pAIC vs pBIC to compare criteria.")
    ),
    mainPanel(width = 9, leafletOutput("map_swipe", height = MAP_H))
  )
)

# ---- Tab 4: Spatial diagnostics (Moran + LISA) -----------------------------
tab_diag <- tabPanel(
  "Spatial diagnostics",
  sidebarLayout(
    sidebarPanel(
      width = 3,
      radioButtons("sd_outcome", "Outcome",
        choices = OUTCOME_CHOICES, selected = "y1", inline = TRUE),
      radioButtons("sd_criterion", "Estimation criterion",
        choices = CRITERION_CHOICES, selected = "pAIC"),
      tags$hr(),
      tags$b("LISA map"),
      selectInput("sd_model", "Model", choices = MODEL_CH_GRP, selected = "0ddd"),
      radioButtons("sd_alpha", "Significance",
        choices = c("p < 0.05" = "0.05", "p < 0.01" = "0.01"),
        selected = "0.05", inline = TRUE),
      checkboxInput("sd_show_legend", "Show legend", TRUE),
      sliderInput("sd_opacity", "Choropleth opacity", 0, 1, 0.85, 0.05),
      helpText("Left: global Moran's I of residuals by model.",
               "Right: local clusters (HH/LL/HL/LH) for the chosen model.")
    ),
    mainPanel(
      width = 9,
      fluidRow(
        column(5, plotOutput("moran_plot", height = "560px")),
        column(7, leafletOutput("map_lisa", height = "560px"))
      )
    )
  )
)

# ---- Tab 5: Model structure (matrices) -------------------------------------
tab_structure <- tabPanel(
  "Model structure",
  sidebarLayout(
    sidebarPanel(
      width = 3,
      helpText(tags$b("Estimated A / R / Λ / Σ."),
               "Conventional models force off-diagonal Λ/Σ to 0 (grey);",
               "the new methodology estimates cross-coupling."),
      radioButtons("ms_criterion", "Estimation criterion",
        choices = CRITERION_CHOICES, selected = "pAIC"),
      selectInput("ms_models", "Models (compare)", choices = MODEL_CH_GRP,
                  selected = c("0ddd", "0111"), multiple = TRUE),
      helpText("Tip: keep 2–3 models for readability.",
               "γ = penalty weight chosen by the criterion (0 = no penalty,",
               "conventional models); d_eff = effective degrees of freedom.")
    ),
    mainPanel(
      width = 9,
      plotOutput("matrix_plot", height = "560px"),
      tags$hr(),
      tags$b(textOutput("fit_title", inline = TRUE)),
      tableOutput("fit_table")
    )
  )
)

ui <- navbarPage(
  title  = "Kansai MGNST viewer — penalized (pAIC/pBIC)",
  header = busy_overlay,
  tab_explore, tab_change, tab_swipe, tab_diag, tab_structure
)

################################################################################
# Server
################################################################################

server <- function(input, output, session) {

  ## ===== Tab 1: Explore ======================================================
  exp_year <- reactive({
    if (input$exp_category %in% RESID_CATS) YEAR_RESID else as.integer(input$exp_year)
  })

  observeEvent(input$exp_category, {
    if (input$exp_category == "residual_diff") return()
    ch <- var_choices(app_data, input$exp_category)
    updateSelectInput(session, "exp_variable", choices = ch, selected = ch[1])
  })

  output$map <- renderLeaflet(base_map())

  # Redraw the choropleth via leafletProxy. We use observeEvent with an explicit
  # trigger list of EVERY input the body reads, so the dependencies are pinned in
  # the trunk. (A bare observe() that read some inputs only inside conditional
  # branches failed to track them reliably here, so those controls did not
  # redraw the map.) Keep input$exp_criterion in this list.
  observeEvent(
    {
      list(input$exp_category, input$exp_variable, input$exp_year,
           input$exp_criterion,
           input$exp_show_choropleth, input$exp_show_legend, input$exp_opacity,
           input$exp_diff_outcome, input$exp_diff_a, input$exp_diff_b)
    },
    {
    req(input$exp_category, input$exp_criterion)
    cat  <- input$exp_category
    crit <- input$exp_criterion

    if (isFALSE(input$exp_show_choropleth)) {
      leafletProxy("map") %>% clearGroup("target") %>% removeControl("legend")
      return()
    }
    yr <- exp_year()

    if (cat == "residual_diff") {
      req(input$exp_diff_outcome, input$exp_diff_a, input$exp_diff_b)
      oc     <- input$exp_diff_outcome
      oc_key <- OUTCOME_KEY[[oc]]
      vals   <- resid_vec(oc_key, input$exp_diff_a, crit) -
                resid_vec(oc_key, input$exp_diff_b, crit)
      req(any(is.finite(vals)))
      var_label <- sprintf("Residual diff: %s (e_%s − e_%s) [%s]",
        OUTCOME_SHORT[[oc]], input$exp_diff_a, input$exp_diff_b, crit)
      sp <- diff_scale(vals)
    } else { # response / covariate / fitted / residual — fixed continuous scale
      req(input$exp_variable)
      if (!(input$exp_variable %in% var_choices(app_data, cat))) return()
      sel <- app_data %>%
        filter(.data$category == cat, .data$variable == input$exp_variable,
               .data$year == yr)
      if (cat %in% RESID_CATS) sel <- sel %>% filter(.data$criterion == crit)
      req(nrow(sel) > 0)
      outcome   <- sel$outcome[1]
      var_label <- if (cat %in% RESID_CATS) sprintf("%s [%s]", sel$label[1], crit)
                   else sel$label[1]
      vals      <- sel$value[match(target$muni_code, sel$muni_code)]
      sp <- scale_for(scales, cat, input$exp_variable, outcome)
    }
    pal <- make_pal(sp)

    labels <- hover_labels(target, vals, var_label)
    popups <- popup_all(app_data, target, yr, crit)

    proxy <- leafletProxy("map") %>%
      clearGroup("target") %>% removeControl("legend") %>%
      addPolygons(data = target, group = "target", layerId = ~muni_code,
        fillColor = pal(vals), fillOpacity = input$exp_opacity,
        color = "#777777", weight = 0.3, opacity = 1, label = labels,
        labelOptions = labelOptions(style = list("font-size" = "12px",
          "max-width" = "210px", "white-space" = "normal")),
        highlightOptions = highlightOptions(weight = 2, color = "#222222",
          fillOpacity = min(1, input$exp_opacity + 0.1), bringToFront = TRUE),
        popup = popups)

    if (isTRUE(input$exp_show_legend)) {
      proxy %>% addLegend(layerId = "legend", position = "bottomright",
        pal = make_legend_pal(sp), values = sp$domain, title = var_label,
        opacity = 0.9, labFormat = labelFormat(digits = 2,
          transform = function(x) sort(x, decreasing = TRUE)))
    }
    },
    ignoreNULL = FALSE
  )

  ## ===== Tab 2: Panel change (Delta) =========================================
  # Whole map built inside renderLeaflet (base + choropleth) so it draws on the
  # first time the tab is shown, with no leafletProxy/hidden-tab timing gap.
  output$map_delta <- renderLeaflet({
    req(input$pc_variable)
    dv  <- delta_vec(input$pc_variable)
    req(any(is.finite(dv$d)))
    lab <- names(PC_CH)[match(input$pc_variable, PC_CH)]
    var_label <- sprintf("Δ %s (%d − %d)", lab, YEAR_RESID, YEAR_BASE)
    sp  <- delta_scale(dv$d); pal <- make_pal(sp)

    labels <- mapply(function(nm, ja, a, b, d)
      HTML(sprintf("<b>%s</b> <span style='color:#888'>(%s)</span><br/>%s<br/>%d: %s &nbsp; %d: %s &nbsp; <b>Δ: %s</b>",
                   nm, ja, var_label, YEAR_BASE, fmt_value(a),
                   YEAR_RESID, fmt_value(b), fmt_value(d))),
      target$muni_name, target$muni_name_ja, dv$v0, dv$v1, dv$d, SIMPLIFY = FALSE, USE.NAMES = FALSE)

    m <- base_map() %>%
      addPolygons(data = target, group = "target",
        fillColor = pal(dv$d), fillOpacity = input$pc_opacity,
        color = "#777777", weight = 0.3, opacity = 1, label = labels,
        highlightOptions = highlightOptions(weight = 2, color = "#222222",
          bringToFront = TRUE))

    if (isTRUE(input$pc_show_legend)) {
      m <- m %>% addLegend(position = "bottomright",
        pal = make_legend_pal(sp), values = sp$domain, title = var_label,
        opacity = 0.9, labFormat = labelFormat(digits = 2,
          transform = function(x) sort(x, decreasing = TRUE)))
    }
    m
  })

  ## ===== Tab 3: Compare (swipe) ==============================================
  output$map_swipe <- renderLeaflet({
    req(input$sw_outcome, input$sw_left, input$sw_right, input$sw_quantity,
        input$sw_left_crit, input$sw_right_crit)
    oc     <- input$sw_outcome
    oc_key <- OUTCOME_KEY[[oc]]
    qn     <- input$sw_quantity
    lv <- vec_for(qn, oc_key, input$sw_left,  input$sw_left_crit)
    rv <- vec_for(qn, oc_key, input$sw_right, input$sw_right_crit)
    req(any(is.finite(lv)), any(is.finite(rv)))

    # Shared scale across the TWO selected maps only (not the global ±rmax,
    # which is dominated by the aspatial models and would wash both halves out).
    if (qn == "residual") {
      m <- max(abs(c(lv, rv)), na.rm = TRUE); if (!is.finite(m) || m == 0) m <- 1
      sp <- list(domain = c(-m, m), palette = "RdBu", diverging = TRUE)
      qlab <- "Residual"
    } else {
      rng <- range(c(lv, rv), na.rm = TRUE)
      if (!all(is.finite(rng)) || diff(rng) == 0) rng <- c(0, 1)
      sp <- list(domain = rng, palette = outcome_pal(oc))
      qlab <- "Fitted"
    }
    pal <- make_pal(sp)

    badge <- function(mid, crit, side) {
      fam <- MODEL_FAMILY[[mid]]
      sprintf("<div style='background:#fff;padding:4px 8px;border-radius:6px;box-shadow:0 1px 4px rgba(0,0,0,.2);font-size:12px'><span class='fam-badge' style='background:%s'>%s</span> %s <span class='crit-badge'>%s</span><br/>Moran I = %.3f</div>",
              FAMILY_COLOR[[fam]], side, MODEL_LABELS[[mid]], crit,
              get_moran(mid, oc, crit))
    }
    mklab <- function(v) lapply(seq_len(nrow(target)), function(i)
      HTML(sprintf("<b>%s</b><br/>%s: %s", target$muni_name[i], qlab, fmt_value(v[i]))))

    # leaflet.extras2 side-by-side clips by PANE and looks the comparison layers
    # up as TILE layers (getLayer("tile", id)). So we add one tile layer per side
    # with a layerId, place each in its own pane, and put each side's polygons in
    # the SAME pane so they are clipped together with that side's tile.
    m <- leaflet(options = leafletOptions(minZoom = 6)) %>%
      addMapPane("left",  zIndex = 410) %>%
      addMapPane("right", zIndex = 420) %>%
      add_basemap(layerId = "leftTile",  pane = "left") %>%
      add_basemap(layerId = "rightTile", pane = "right") %>%
      addPolygons(data = target, group = "left",
        options = pathOptions(pane = "left"),
        fillColor = pal(lv), fillOpacity = input$sw_opacity,
        color = "#777777", weight = 0.3, label = mklab(lv)) %>%
      addPolygons(data = target, group = "right",
        options = pathOptions(pane = "right"),
        fillColor = pal(rv), fillOpacity = input$sw_opacity,
        color = "#777777", weight = 0.3, label = mklab(rv)) %>%
      addPolylines(data = kansai_borders, color = "#888888",
                   weight = 0.8, opacity = 0.6) %>%
      addSidebyside(layerId = "cmp", leftId = "leftTile", rightId = "rightTile") %>%
      addControl(HTML(badge(input$sw_left,  input$sw_left_crit,  "LEFT")),
                 position = "topleft") %>%
      addControl(HTML(badge(input$sw_right, input$sw_right_crit, "RIGHT")),
                 position = "topright") %>%
      fitBounds(bb[1], bb[2], bb[3], bb[4])

    if (isTRUE(input$sw_show_legend)) {
      m <- m %>% addLegend(position = "bottomright", pal = make_legend_pal(sp),
        values = sp$domain, title = sprintf("%s — %s", qlab, OUTCOME_SHORT[[oc]]),
        opacity = 0.9, labFormat = labelFormat(digits = 2,
          transform = function(x) sort(x, decreasing = TRUE)))
    }
    m
  })

  ## ===== Tab 4: Spatial diagnostics ==========================================
  output$moran_plot <- renderPlot({
    req(input$sd_outcome, input$sd_criterion)
    build_moran_plot(moran_df, input$sd_outcome, input$sd_criterion)
  })

  # Whole map built inside renderLeaflet (see Panel change note above).
  output$map_lisa <- renderLeaflet({
    req(input$sd_model, input$sd_outcome, input$sd_alpha, input$sd_criterion)
    alpha <- as.numeric(input$sd_alpha)
    d <- lisa_df %>% filter(.data$model == input$sd_model,
                            .data$outcome == input$sd_outcome,
                            .data$criterion == input$sd_criterion)
    req(nrow(d) > 0)
    d$cluster <- lisa_cluster(d$quadrant, d$p_value, alpha)
    cl <- d$cluster[match(target$muni_code, d$muni_code)]
    pv <- d$p_value[match(target$muni_code, d$muni_code)]
    pal_fac <- colorFactor(LISA_COLORS, levels = LISA_LEVELS, na.color = NA_COLOR)

    labels <- mapply(function(nm, ja, c, p)
      HTML(sprintf("<b>%s</b> <span style='color:#888'>(%s)</span><br/>LISA: %s<br/>p = %s",
                   nm, ja, LISA_LABELS[as.character(c)], fmt_value(p))),
      target$muni_name, target$muni_name_ja, cl, pv, SIMPLIFY = FALSE, USE.NAMES = FALSE)

    m <- base_map() %>%
      addPolygons(data = target, group = "target",
        fillColor = pal_fac(as.character(cl)), fillOpacity = input$sd_opacity,
        color = "#777777", weight = 0.3, opacity = 1, label = labels,
        highlightOptions = highlightOptions(weight = 2, color = "#222222",
          bringToFront = TRUE))

    if (isTRUE(input$sd_show_legend)) {
      m <- m %>% addLegend(position = "bottomright",
        colors = unname(LISA_COLORS), labels = unname(LISA_LABELS),
        title = sprintf("LISA cluster — %s, %s", input$sd_model,
                        input$sd_criterion),
        opacity = 0.9)
    }
    m
  })

  ## ===== Tab 5: Model structure ==============================================
  output$matrix_plot <- renderPlot({
    req(input$ms_models, input$ms_criterion)
    build_matrix_plot(coef_df, input$ms_models, input$ms_criterion)
  })

  output$fit_title <- renderText({
    req(input$ms_criterion)
    sprintf("Model fit (sorted by %s)", input$ms_criterion)
  })

  output$fit_table <- renderTable({
    req(input$ms_criterion)
    crit <- input$ms_criterion
    d <- fit_df
    d$crit_val <- d[[crit]]
    d <- d[order(d$crit_val), ]
    out <- data.frame(
      Model    = unname(MODEL_LABELS[d$model]),
      Family   = unname(FAMILY_LABEL[d$family]),
      `n par`  = as.integer(d$n_params),
      d_eff    = d[[paste0("d_eff_", crit)]],
      `γ`      = formatC(d[[paste0("gamma_", crit)]], format = "g", digits = 4),
      logLik   = d[[paste0("loglik_", crit)]],
      crit     = d$crit_val,
      delta    = d[[paste0("delta_", crit)]],
      `R²`     = d$R2,
      check.names = FALSE, stringsAsFactors = FALSE
    )
    names(out)[names(out) == "crit"]  <- crit
    names(out)[names(out) == "delta"] <- paste0("Δ", crit)
    out
  }, digits = 2, striped = TRUE, width = "100%")
}

shinyApp(ui, server)
