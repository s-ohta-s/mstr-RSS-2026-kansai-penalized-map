# Kansai MGNST viewer — penalized (pAIC/pBIC)

Interactive Shiny + leaflet viewer of the **penalized** (pAIC / pBIC) estimation results for
the 198 Kansai municipalities. Theme: *Visual Comparison of New Spatial Methodology with
Conventional Individual Regression Models.*

The app **only reads pre-computed `.rds` files**; it performs no estimation and no spatial
statistics at runtime. Everything is built by the scripts in `prep/`.

**Live demo:** <https://01a0ccab-923c-ef0c-b3be-aa4fc076178b.share.connect.posit.cloud/>

## What changed from the non-penalized viewer

The layout (five tabs) follows the earlier non-penalized viewer
([s-ohta-s/spatial-interactions](https://github.com/s-ohta-s/spatial-interactions)). Changes:

- **Criterion axis is back.** Each penalized model has its own optimal γ for pAIC and for
  pBIC, so fitted values and residuals can differ between the two. Every tab has a criterion
  selector; the swipe tab picks a criterion per side, so the same model can be compared under
  pAIC and pBIC.
- **Innovation residuals.** Residual maps, Moran's I and LISA use the innovation residuals
  `y*_resid` of the estimation. The reduced-form residuals `y*_resid_reduced` are not used.
- **New dataset and years** (*Data_I*): responses are the share of population aged 65+ and the
  share of secondary-industry workers (both standardized); four new covariates; years FY2015
  and FY2020 (fitted values and residuals: FY2020 only).
- **Estimation results** come from the `real_v2` run (R 4.6.1).
- **Basemap key is optional.** With `CARTO_BASEMAP_KEY` the keyed CARTO URL is used; without it
  the app falls back to `providers$CartoDB.Positron` instead of stopping (see below).
- **Constants in one place.** `R/constants.R` holds the variables, labels, years, criteria and
  model tables; both `prep/` and the app source it.

## Tabs

1. **Explore** — one map of a response, covariate, fitted value, residual, or residual
   difference `e_A − e_B` (any two models, same criterion). Hover for the value, click for a
   table of all variables and the per-model residuals / fitted values of the selected
   criterion.
2. **Panel change** — Δ = value(2020) − value(2015) for responses and covariates (diverging,
   0-centered).
3. **Compare (swipe)** — drag the handle to wipe between two maps of residuals or fitted
   values. Each side has its own model and criterion (default: `0ddd` pAIC vs `0111` pAIC).
   Both sides share one color scale taken from the two maps. The badges show the family,
   criterion and residual Moran's I.
4. **Spatial diagnostics** — global Moran's I of the residuals by model (new vs conventional)
   and a LISA cluster map (HH / LL / HL / LH at p < 0.05 or 0.01), for the selected criterion.
5. **Model structure** — estimated A / R / Λ / Σ matrices of the selected models, and a model
   fit table sorted by the selected criterion (with d_eff, γ, log-likelihood and Δ).
   Significance marks use Ψ-based standard errors; where those are not available (R, Λ) the
   Hessian-based ones are used. Σ and the conventional models' R / Λ have no standard error
   and no mark.

## Color scales

Matching `reference/create_heatmaps_en.r`:

- y1 (aged 65+) = viridis *plasma*, y2 (secondary industry) = viridis *mako*, covariates =
  *viridis*; residuals and differences = *RdBu*, 0-centered (low = red, high = blue).
- Responses and fitted values share one range per outcome, taken over the observed values and
  the fitted values of all 11 models × 2 criteria. Residuals share one symmetric range over all
  models, outcomes and criteria. These ranges are set by the worst-fitting model (`000d`), so
  the residuals of the well-fitting models appear in pale colors in **Explore**. The swipe tab
  uses the range of its two maps only.

## Run locally

From the repository root:

```r
shiny::runApp("shiny_app")
```

Runtime packages: `shiny`, `leaflet`, `leaflet.extras2`, `sf`, `dplyr`, `htmltools`,
`ggplot2`, `RColorBrewer`, `viridisLite`. The prep scripts also need `spdep`.

### Basemap key

Set `CARTO_BASEMAP_KEY` in your user `.Renviron` (on Windows: `Documents\.Renviron`), or as a
secret on Posit Connect Cloud. Never commit it. At startup the app prints which basemap it uses
(never the key):

```
[basemap] CARTO Positron (API key from CARTO_BASEMAP_KEY)
[basemap] CARTO_BASEMAP_KEY not set -> fallback providers$CartoDB.Positron
```

Without a key, CARTO serves tiles with an **"API KEY REQUIRED" watermark**; everything else
works.

## Files

```
shiny_app/
├── app.R                      # UI (navbarPage, 5 tabs) + server
├── R/
│   ├── constants.R            # variables, years, criteria, model tables (shared with prep/)
│   └── helpers.R              # scales, labels, popups, plot builders
├── prep/
│   ├── build_app_data.R       # data/real + output/real -> app_data.rds
│   ├── build_spatial_diag.R   # W + residuals -> spatial_diag.rds (Moran's I, LISA)
│   └── build_matrices.R       # comparison tables -> coef_matrices.rds
├── data/
│   ├── muni_sf.rds            # polygons: 198 target + 7 background, EPSG:4326 (reused)
│   ├── pref_borders.rds       # 13 prefecture outlines (reused)
│   ├── app_data.rds           # long table: response / covariate / fitted / residual
│   ├── spatial_diag.rds       # list(moran, lisa)
│   └── coef_matrices.rds      # list(coef, fit, gamma)
├── DESIGN_penalized.md        # design notes (Japanese)
├── .renvignore                # keeps prep/ out of renv / deployment
└── renv.lock
```

### `.rds` contents

| File | Rows | Key columns |
|---|---|---|
| `app_data.rds` | 19,800 | `muni_code, year, category, variable, outcome, model, criterion, value, label`. `criterion` is `pAIC` / `pBIC` for fitted and residual rows, `NA` for responses and covariates |
| `spatial_diag.rds$moran` | 44 | `model, outcome, criterion, morans_i, variance, z, p_analytic, p_value` (permutation, 999 sims, seed 20260614) |
| `spatial_diag.rds$lisa` | 8,712 | `muni_code, model, outcome, criterion, local_i, quadrant, p_value` (NS is applied in the app) |
| `coef_matrices.rds$coef` | 352 | `model, criterion, matrix, row, col, estimate, se, signif, signif_source, present` |
| `coef_matrices.rds$fit` | 11 | `n_params, d_eff_*, loglik_*, AIC, pAIC, BIC, pBIC, gamma_*, R2, delta_*` |
| `coef_matrices.rds$gamma` | 6 | optimal γ, d_eff and log-likelihood per criterion |

## Rebuilding the data

From the repository root (see the root [README](../README.md)):

```sh
Rscript shiny_app/prep/build_app_data.R
Rscript shiny_app/prep/build_spatial_diag.R
Rscript shiny_app/prep/build_matrices.R
```

Inputs: `data/real/transformed_data198_for_R.csv`, `region_id_mapping198.csv`
(columns `region, code, name`), `W_row_standardized_198_for_R.csv`, and
`output/real/{comparison_table_p*.csv, model_selection_summary.csv,
gamma_optimal_summary.csv, fitted_residuals_p*/}`. Set `MSNR_REPO_ROOT` to run from elsewhere.

## Deployment (Posit Connect Cloud)

The app is deployed with Posit Publisher from this folder, using
`.posit/publish/Kansai MGNST viewer — penalize-O1E2.toml` (a separate deployment from the
non-penalized viewer). To update it, open `shiny_app/` in Positron and click **Deploy** again.

- `files` must list `app.R`, `R`, `data`, `renv.lock` and `.renvignore` (`prep/` is not needed
  at runtime). Publisher's default list contained only `app.R`, which would fail at startup.
- For Connect Cloud, `[r]` accepts only `version`; `package_file` / `package_manager` are
  rejected by the schema. `renv.lock` is uploaded as a file.
- `CARTO_BASEMAP_KEY` is a secret (only its name is in the `.toml`).
- `.posit/publish/deployments/` (per-machine deployment records) is git-ignored.
