# mstr-RSS-2026-kansai-penalized-map

**Interactive maps of penalized (pAIC / pBIC) multivariate spatio-temporal regression results
for the 198 municipalities of the Kansai region, Japan.**

This repository ships a Shiny viewer for the estimation results of the multivariate
spatio-temporal models with **penalized model selection** (pAIC / pBIC). It compares the new
multivariate spatial methodology (MGNS family) with conventional, individually fitted
regression models, on maps and diagnostics. It is the penalized counterpart of the earlier
non-penalized viewer ([s-ohta-s/spatial-interactions](https://github.com/s-ohta-s/spatial-interactions)).

The app only **reads pre-computed files**; it does not run any estimation.

---

## Overview

- **Area**: all 198 municipalities of Shiga, Kyoto, Osaka, Hyogo, Nara and Wakayama.
- **Periods**: FY2015 (`time = 1`) and FY2020 (`time = 2`). Fitted values and residuals exist
  for FY2020 only.
- **Variables** (dataset *Data_I*; definitions in
  [`data/real/data_source_and_processing_en.md`](data/real/data_source_and_processing_en.md)):

  | Column | Role | Label in the app |
  |---|---|---|
  | `y1` | Response 1 | Population Aged 65+ (%, standardized) |
  | `y2` | Response 2 | Secondary-Industry Workers (%, standardized) |
  | `x_common1` | Covariate (both equations) | Log Taxable Income per Taxpayer |
  | `x_common2` | Covariate (both equations) | Log(1 + Foreign Population per 100k) |
  | `x_specific1_1` | Covariate (y1 equation) | Net In-migration Rate, Japanese migrants (%) |
  | `x_specific2_1` | Covariate (y2 equation) | Commuters from Other Municipalities (%) |

- **Models** (11):

  | Model ID | Name | Family |
  |---|---|---|
  | `1111` | MGNS | New methodology |
  | `0111` | MSEM | New methodology |
  | `1011` | MSAR | New methodology |
  | `1101` | MGNS | New methodology |
  | `0101` | MSEM | New methodology |
  | `1001` | MSAR | New methodology |
  | `d0dd` | IndSAR | Conventional (individual) |
  | `0ddd` | IndSEM | Conventional (individual) |
  | `dddd` | IndGNS | Conventional (individual) |
  | `0011` | VARX | Conventional (individual) |
  | `000d` | IndReg | Conventional (individual) |

- **Penalized selection**: each of the six new-methodology models has its own optimal penalty
  weight γ for each criterion, so their fitted values and residuals can differ between pAIC
  and pBIC. In this estimation they differ for `1011`, `1101`, `0101`, `1001` and are identical
  for `1111`, `0111` (same γ in both criteria). Conventional models are not penalized (γ = 0).
- **Residuals** are the **innovation residuals** of the estimation (`y*_resid`). The
  reduced-form residuals (`y*_resid_reduced`) are not used.

---

## Repository structure

```
.
├── README.md                 # this file
├── LICENSE                   # MIT (code)
├── PLAN.md                   # work plan and work log (Japanese)
├── data/real/                # input data (Data_I) and its documentation
│   ├── transformed_data198_for_R.csv
│   ├── region_id_mapping198.csv          # region (1..198) <-> municipal code, name
│   ├── W_row_standardized_198_for_R.csv  # row-standardized spatial weights
│   ├── Data_I_ver2.txt
│   ├── DATA_SOURCE_AND_PROCESSING.md     # data definition (Japanese)
│   └── data_source_and_processing_en.md  # data definition (English)
├── output/real/              # estimation results (penalized, pAIC / pBIC)
│   ├── comparison_table_pAIC.csv / comparison_table_pBIC.csv
│   ├── model_selection_summary.csv
│   ├── gamma_optimal_summary.csv
│   └── fitted_residuals_pAIC/ , fitted_residuals_pBIC/   # 11 CSVs each
├── reference/                # static-map script and its console log (used for number checks)
└── shiny_app/                # the viewer (see shiny_app/README.md)
    ├── app.R
    ├── R/constants.R         # variables, years, criteria, model tables (shared by prep and app)
    ├── R/helpers.R
    ├── prep/                 # build the .rds files below from data/ and output/
    └── data/                 # pre-computed .rds read by the app
```

---

## Quick start

Requirements: R (tested with 4.6.1) and the packages below.

```r
install.packages(c("shiny", "leaflet", "leaflet.extras2", "sf", "dplyr",
                   "htmltools", "ggplot2", "RColorBrewer", "viridisLite"))
# prep only (not needed to run the app):
install.packages("spdep")
```

From the repository root:

```r
shiny::runApp("shiny_app")
```

The `.rds` files in `shiny_app/data/` are committed, so the app runs without the prep step.

### Basemap key (optional)

The basemap is CARTO Positron. With an API key the app requests keyed CARTO tiles; without one
it falls back to leaflet's `CartoDB.Positron` provider. **Without a key, CARTO returns tiles
with an "API KEY REQUIRED" watermark**; the app itself works normally.

To use a key, put it in your user `.Renviron` (never commit it; `.Renviron` is in
`.gitignore`):

```
CARTO_BASEMAP_KEY=<your key>
```

On Windows, R's user `.Renviron` lives in `Documents\` (R's `HOME`). On Posit Connect Cloud,
set `CARTO_BASEMAP_KEY` as a secret.

---

## Reproducing the app data

The three prep scripts rebuild the app's data from `data/real/` and `output/real/`. Run them
from the repository root:

```sh
Rscript shiny_app/prep/build_app_data.R      # -> shiny_app/data/app_data.rds
Rscript shiny_app/prep/build_spatial_diag.R  # -> shiny_app/data/spatial_diag.rds
Rscript shiny_app/prep/build_matrices.R      # -> shiny_app/data/coef_matrices.rds
```

Each script prints sanity checks and stops if one fails. `build_spatial_diag.R` checks the
residual Moran's I of models `1111` and `dddd` (pAIC) against `reference/`.

The municipality polygons (`muni_sf.rds`, `pref_borders.rds`) do not depend on the estimation
and are reused from the earlier viewer; they are not rebuilt here.

---

## Data sources

- **Socioeconomic statistics** (the variables above) — processed from official statistics:

  > Created from the Statistics Bureau of Japan, "System of Social and Demographic
  > Statistics: Municipal Data, Social Indicators (Adjusted for Municipal Mergers and Boundary
  > Changes)," A: Population and Households (Table ID 0000020301), D: Administrative Base
  > (Table ID 0000020304), and F: Labor (Table ID 0000020306), Portal Site of Official
  > Statistics of Japan (e-Stat), retrieved July 3, 2026.

  See [`data/real/data_source_and_processing_en.md`](data/real/data_source_and_processing_en.md)
  for the indicator codes and transformations. When you reuse the processed data, cite the
  original data above and this repository (with its commit ID and access date), and follow the
  [e-Stat Terms of Use](https://www.e-stat.go.jp/terms-of-use).
- **Municipality polygons** — National Land Numerical Information, Administrative Areas (N03),
  Ministry of Land, Infrastructure, Transport and Tourism:
  <https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N03-v3_1.html>. The bundled polygons are
  simplified and clipped to Kansai and its neighbouring prefectures. (The N03 release used to
  build them is not recorded in this repository; the reference static maps use
  `N03-20210101`.)
- **Basemap** — © [OpenStreetMap](https://www.openstreetmap.org/copyright) contributors,
  © [CARTO](https://carto.com/attributions) (shown on the map).

---

## License

- **Code** (`shiny_app/`, including `prep/`) — [MIT License](LICENSE).
- **Data** (`data/`, `output/`, `shiny_app/data/`) — processed from the sources above. Their
  reuse is governed by the original providers' terms (e-Stat Terms of Use for the statistics,
  the National Land Numerical Information terms for the polygons). Cite the sources as shown
  above.

## Author

[@s-ohta-s](https://github.com/s-ohta-s)
