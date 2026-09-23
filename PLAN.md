# 罰則あり推定 関西 198 市町村 Shiny マップ — 作成計画（Opus 5.5 引き継ぎ用・確定版）

作成: 2026-09-23 Claude Fable 5.1（調査・計画のみ。実装は本計画に基づき別セッションで行う）
確定: 2026-09-23 ユーザー回答（§1）を反映

本フォルダ `C:\Users\s.ota\Desktop\残差地図_罰則有\` を作業ルートとし、
罰則あり（pAIC / pBIC）推定の結果を閲覧する Shiny 地図アプリを、
前回の罰則なし版アプリ（`20260613_noPenalty/shiny_app/`）を土台に作る。

---

## 0. 調査で確定した事実（実装者はここを前提にしてよい）

### 0.1 参照するソースフォルダ（すべて **読み取り専用**。編集・上書き禁止）

| フォルダ | 役割 | 使うもの |
|---|---|---|
| `G:\マイドライブ\01study\04_空間統計\20260613_noPenalty\` | 前回（罰則なし）の Shiny アプリ一式 | `shiny_app/`（app.R, R/helpers.R, prep/*.R, data/*.rds, README.md, DESIGN*.md, .renvignore, renv.lock）、`.gitignore` |
| `G:\マイドライブ\01study\04_空間統計\20260705_Kansai_Estimation\` | 罰則あり推定の入力データと推定結果 | `data/real/`、**`real_v2/`**（最新の推定出力） |
| `G:\マイドライブ\01study\04_空間統計\20260709_Kansai_Map\` | 静止画（ggplot）版マップ | `create_heatmaps_en.r`（配色の一次参照）、`heatmaps_Data_I_v2/console.txt`（数値照合用） |
| `G:\マイドライブ\01study\04_空間統計\20260702_RealDataKansai\Data_I\` | データ定義書 | `DATA_SOURCE_AND_PROCESSING.md`（日本語）、`data_source_and_processing_en.md`（英語） |

### 0.2 推定出力は `real_v2/` が正（確定）

| 場所 | 日付 | 残差 CSV の列 | 備考 |
|---|---|---|---|
| `output/real/` | 2026-07-19 | 8 列（`region,time,y1_obs,y1_pred,y1_resid,y2_obs,y2_pred,y2_resid`） | 旧形式。使わない |
| **`real_v2/`** | **2026-07-24** | **10 列**（上記 + `y1_resid_reduced`, `y2_resid_reduced`） | R 4.6.1 で再実行。`20260709_Kansai_Map/fitted_residuals_p*` と同一（cmp 確認済） |

補足（旧形式について）: 旧形式の `resid` 列も新形式の `resid`（イノベーション残差）と**数値は完全一致**する。
違いは `pred` の定義で、旧 = `obs − ε̂`（恒等式で作った値）、新 = 誘導形トレンド予測 ŷ。
したがって旧形式を混ぜると**推定値地図だけ**がずれる。いずれにせよ real_v2 のみ使う。

### 0.3 残差 CSV の列の意味（`R/inference_and_export.r` §5 のコメントより）

- `pred` = 誘導形トレンド予測 ŷ（d = 0）。OLS/VARX/MSEM: ŷ = Xβ̂、MSAR/MGNS: ŷ = (I − R̂⊗W)⁻¹Xβ̂
- `resid` = **イノベーション残差 ε̂**（尤度・Σ 推定の残差）。**アプリで地図化・Moran's I・LISA に使うのはこれだけ（確定）**
- `resid_reduced` = `obs − pred`（誘導形残差）。**アプリでは使わない**（読み込んでも捨てる）
- 検証済: `obs − pred == resid_reduced`（1111）、非空間モデル 000d では `resid == resid_reduced`

### 0.4 pAIC / pBIC で fitted / residual が異なるモデル

| モデル | pAIC と pBIC の残差ファイル | 理由 |
|---|---|---|
| 1011, 1101, 0101, 1001 | **異なる**（最大差 0.06〜0.35） | 最適 γ が基準で異なる（`gamma_optimal_summary.csv`） |
| 1111, 0111 | 同一 | γ が両基準で同じ |
| d0dd, 0ddd, dddd, 0011, 000d | 同一 | 従来モデルは罰則なし（γ = 0） |

→ 前回アプリで撤去した **criterion（pAIC / pBIC）軸を復活**させる。

### 0.5 データセット Data_I（前回と別の変数・年）— 定義書で確定

`transformed_data198_for_R.csv`（`20260705_Kansai_Estimation/data/real/`）。列名は前回と同じだが中身は別変数。
定義書 `data_source_and_processing_en.md` §3.1 に基づく英語ラベル:

| 列 | key（app 内部名） | 英語ラベル（UI 表示） | 備考 |
|---|---|---|---|
| `y1` | `aging65` | Population Aged 65+ (%, standardized) | e-Stat #A03506、396 観測プールで標準化 |
| `y2` | `secondary_ind` | Secondary-Industry Workers (%, standardized) | e-Stat #F01202 |
| `x_common1` | `ln_taxable_income` | Log Taxable Income per Taxpayer | ln(千円) |
| `x_common2` | `ln_foreign_pop` | Log(1 + Foreign Population per 100k) | ln(x+1) |
| `x_specific1_1` | `net_inmigration` | Net In-migration Rate, Japanese migrants (%) | 無変換。y1 式のみ |
| `x_specific2_1` | `commuter_inflow` | Commuters from Other Municipalities (%) | 無変換。y2 式のみ。最大 231 で分布が偏る |
| `time == 1 / 2` | — | **FY2015 / FY2020** | UI では "2015" / "2020" 表記でよい |

定義書の記述統計（396 観測、SD は母標準偏差）はゲートに使える:
y1 [-1.8887, 3.5265]、y2 [-2.1456, 2.9951]、x_common1 [7.7212, 8.7828]、x_common2 [0, 8.6415]、x_specific1_1 [-5.79, 1.41]、x_specific2_1 [7.8, 231.1]。

### 0.6 地理データ・W・mapping は前回のものをそのまま流用できる

- `20260613_noPenalty/shiny_app/data/muni_sf.rds`（198 対象 + 背景 7 県、EPSG:4326、列 `muni_code, muni_name, muni_name_ja, is_target, pref_code, pref_name`）と `pref_borders.rds`（13 府県境）は推定に依存しない
- `W_row_standardized_198_for_R.csv` と `region_id_mapping198.csv` は **前回プロジェクトのものとバイト単位で同一**（cmp 確認済）。mapping の `code` 先頭 5 桁は `muni_sf` の対象 198 コードと全件一致
- → シェープファイル（N03）や `rmapshaper` は不要

### 0.7 マッピングファイルの列名の罠

| ファイル | 列名 |
|---|---|
| `20260705_Kansai_Estimation/data/real/region_id_mapping198.csv` | `region, code, name`（前回 prep が期待する形式。**これを使う**） |
| `20260709_Kansai_Map/region_id_mapping198_en.csv` | `region_id, code_6, municipality`（名前は日本語のまま。使わない） |

英語市町村名は `muni_sf.rds` 側にしかない。

### 0.8 推定サマリ CSV の構造変化（build_matrices.R に影響）

- `comparison_table_pAIC.csv` / `comparison_table_pBIC.csv` の 2 本
  - 行: β 各種、`A[i,j]`, `R[i,j]`, `Lambda[i,j]`, `Sigma[1,1]/[1,2]/[2,2]`, `loglik, AIC, BIC, pAIC, pBIC, d_eff, gamma, n_params, pseudo_R2, adj_R2, pseudo_R2_y1, pseudo_R2_y2`
  - 列: モデルごとに `estimate.<ID>, se_psi.<ID>, se_hessian.<ID>, signif_psi.<ID>, signif_hessian.<ID>`
  - **R / Lambda は `se_psi` が空で `se_hessian` のみ** → `signif_psi` が空なら `signif_hessian` を使う（確定、注記付き）
- `model_selection_summary.csv`: `model_id, n_params, d_eff_pAIC, d_eff_pBIC, loglik_pAIC, loglik_pBIC, AIC, pAIC, BIC, pBIC, gamma_pAIC, gamma_pBIC, R2_mean, adj_R2_mean, delta_AIC, delta_pAIC, delta_BIC, delta_pBIC`。**`model_type` 列が無い** → 定数表から補う
- `gamma_optimal_summary.csv`（6 罰則モデルの最適 γ と d_eff）

### 0.9 R 環境（このマシン）

| 項目 | 状態 |
|---|---|
| R 4.6.1 | `C:\Program Files\R\R-4.6.1\bin\Rscript.exe`。shiny 1.14, leaflet 2.2.3, leaflet.extras2 1.3.2, sf 1.1.1, dplyr 1.2.1, ggplot2 4.0.3, spdep 1.4.2, viridisLite, RColorBrewer, renv 1.2.4 すべて有り |
| R 4.3.3 | `C:\R-4.3.3\`（旧。前回の launch.json はこちら → 4.6.1 に差し替え） |
| PATH | `Rscript` は PATH に無い。**必ずフルパスで呼ぶ** |
| 注意 | G: ドライブ上で sf + readRDS を `Rscript -e` の複数行で実行したら一度セグフォルトした。**作業は C: の本フォルダで、R コードはファイルにして `Rscript file.R`** |

### 0.10 CARTO ベースマップと API キーの扱い（確定）

- ユーザーは `CARTO_BASEMAP_KEY` を持っている。**キーはリポジトリに絶対に入れない**
- ローカル: ユーザー自身が `%USERPROFILE%\Documents\.Renviron`（R の HOME。現状は存在しない）に `CARTO_BASEMAP_KEY=...` を書く。実装者はキーを扱わず、「設定してください」と依頼するだけ
- `.gitignore` に `.Renviron` を追加する
- アプリは「キーがあれば鍵付き URL、無ければ `providers$CartoDB.Positron` にフォールバック」とする（クローンした第三者も動かせるように）。前回の `stop()` は撤去
- Posit Connect Cloud へ出す場合は secret として設定（前回と同じ）

---

## 1. 確認事項と回答（2026-09-23 確定）

| # | 論点 | 決定 |
|---|---|---|
| Q1 | 正とする推定出力 | **`real_v2/`** |
| Q2 | 地図化する残差 | **イノベーション残差 `resid` のみ**。`resid_reduced` は UI に出さない（residual-type 軸は作らない） |
| Q3 | 年・ラベル | **FY2015 / FY2020、§0.5 のラベル**。定義書は `20260702_RealDataKansai/Data_I/` の 2 本 |
| Q4 | CARTO キー | **持っている**。§0.10 の運用 |
| Q5 | git | **新規リポジトリを作り、本フォルダをその作業ツリーにする**。手順は §2.1。リポジトリ名は着手前に決める（提案は §2.1） |
| Q6 | 有意性 | **Ψ-SE が無い場合は Hessian-SE の有意性を使い、注記を付ける** |
| Q7 | criterion 切替 | **各タブ共通セレクタ + Swipe のみ左右独立**（同一モデルの pAIC vs pBIC 比較を可能に） |
| Q8 | UI 言語・タイトル | **英語。タイトルに "penalized (pAIC/pBIC)"** |

---

## 2. 成果物のフォルダ構成（本フォルダ内に自己完結させる）

```
残差地図_罰則有/                 ← git の作業ツリー（§2.1）
├── PLAN.md                      # 本計画（作業ログを末尾に追記）
├── README.md                    # プロジェクト概要・再現手順・データ出典（Phase 5）
├── LICENSE                      # 前回と同じ MIT（コード）。データは README で出典・条件を明記
├── .gitignore                   # 前回のものを流用 + `.Renviron`
├── .claude/launch.json          # R 4.6.1 で shiny_app を起動（git 管理外）
├── data/real/                   # 入力（20260705_Kansai_Estimation/data/real/ からコピー）
│   ├── transformed_data198_for_R.csv
│   ├── region_id_mapping198.csv          # 列: region, code, name
│   ├── W_row_standardized_198_for_R.csv
│   ├── Data_I_ver2.txt
│   ├── DATA_SOURCE_AND_PROCESSING.md     # 20260702_RealDataKansai/Data_I/ から
│   └── data_source_and_processing_en.md  # 同上
├── output/real/                 # 推定結果（20260705_Kansai_Estimation/real_v2/ からコピー）
│   ├── comparison_table_pAIC.csv
│   ├── comparison_table_pBIC.csv
│   ├── model_selection_summary.csv
│   ├── gamma_optimal_summary.csv
│   ├── fitted_residuals_pAIC/fitted_residuals_<ID>.csv   # 11 本
│   └── fitted_residuals_pBIC/fitted_residuals_<ID>.csv   # 11 本
├── reference/                   # 照合用（読むだけ）
│   ├── create_heatmaps_en.r                # 20260709_Kansai_Map/ から
│   └── heatmaps_Data_I_v2_console.txt      # 20260709_Kansai_Map/heatmaps_Data_I_v2/console.txt
└── shiny_app/
    ├── app.R
    ├── R/
    │   ├── constants.R          # ★新規: 変数定義・年・モデル表・criterion を一元化（prep と app が共に source）
    │   └── helpers.R
    ├── prep/
    │   ├── build_app_data.R
    │   ├── build_spatial_diag.R
    │   └── build_matrices.R
    ├── data/
    │   ├── muni_sf.rds          # 前回から流用（再生成しない）
    │   ├── pref_borders.rds     # 前回から流用（再生成しない）
    │   ├── app_data.rds         # prep で生成
    │   ├── spatial_diag.rds     # prep で生成
    │   └── coef_matrices.rds    # prep で生成
    ├── README.md
    ├── DESIGN_penalized.md      # 前回 DESIGN.md / DESIGN_beta.md からの差分仕様
    ├── .renvignore              # 前回から流用（prep/ を除外）
    └── renv.lock                # 前回のものを初期値としてコピー → Phase 6 で snapshot し直す
```

**コピーしないもの**: `shiny_app/.posit/`（前回アプリの配布 ID を含む。コピーすると前回の公開アプリを上書きしうる）、`shiny_app_alpha/`、`old_app.R`、`DESIGN.md`/`DESIGN_beta.md` 本体（参照のみ）。

### 2.1 git リポジトリの作り方（Q5）

「新規リポジトリを本フォルダにクローンする」形にしたい場合、空リポジトリは既存の非空フォルダへ `git clone` できないので、実際の手順は次のどちらか。

- **A（推奨）: 本フォルダで `git init` → GitHub に空リポジトリを作成 → `git remote add origin` → push。** フォルダ名（日本語）とリポジトリ名が違っても問題ない
- B: 先に GitHub でリポジトリを作り、Desktop に英語名でクローンし、そこを作業フォルダにする（本計画のパスを読み替える）

リポジトリ名は README・Posit のタイトル・引用表記に出るので **Phase 0 で決める**。前回が `s-ohta-s/spatial-interactions` なので、候補は `spatial-interactions-penalized`（前回との対応が明確）または `kansai-penalized-map`。
公開設定は最初は **private** で作り、公開時に README のデータ出典（定義書 §8 の引用例、N03 の出典）が揃っていることを確認する。

---

## 3. 作業手順

各 Phase の末尾に **ゲート（確認項目）** を置く。ゲートを通過してから次へ進む。

### Phase 0. 環境準備・リポジトリ（30 分）

1. `.claude/launch.json` を作成（前回のものを流用し Rscript を 4.6.1 に変更、ポート 7897）:
   ```json
   {
     "version": "0.0.1",
     "configurations": [
       {
         "name": "shiny-penalized",
         "runtimeExecutable": "C:\\Program Files\\R\\R-4.6.1\\bin\\Rscript.exe",
         "runtimeArgs": ["-e", "shiny::runApp('shiny_app', port=7897, host='127.0.0.1', launch.browser=FALSE)"],
         "port": 7897
       }
     ]
   }
   ```
2. `.gitignore` を前回から流用（`shiny_app_alpha/` の行は削除、`.Renviron` を追加）。`git init`、リポジトリ名を決めて remote 設定（§2.1）。
3. ユーザーに `Documents\.Renviron` へ `CARTO_BASEMAP_KEY` を設定してもらう（実装者はキーを見ない）。
4. **ゲート**: `"C:\Program Files\R\R-4.6.1\bin\Rscript.exe" -e "packageVersion('shiny'); Sys.getenv('CARTO_BASEMAP_KEY') != ''"` が `TRUE` を返す。

### Phase 1. データ配置（15 分）

1. §2 のとおりコピー（PowerShell `Copy-Item`）。`output/real/` は **`real_v2/` から**。
2. **ゲート（必ず数で確認）**:
   - `fitted_residuals_pAIC/` と `_pBIC/` に各 11 ファイル、各 198 行 + ヘッダ、列数 10、`time` はすべて 2
   - `region_id_mapping198.csv` の列名が `region,code,name`、198 行
   - `transformed_data198_for_R.csv` が 396 行、各列の min/max が §0.5 の定義書の値と一致
   - コピー元との `Get-FileHash` 一致
   - `muni_sf.rds` を読んで `is_target` TRUE 198 / FALSE 7、対象コードが mapping の `substr(code,1,5)` と全件一致
3. 最初のコミット（データ + 流用 rds + PLAN.md）。

### Phase 2. prep スクリプトの移植・改修（半日）

共通方針: **アプリは読むだけ**。推定・空間統計はすべて prep で事前計算して `.rds` に落とす。
パスは `REPO_ROOT`（既定 `getwd()` = 本フォルダ）基準。本フォルダで `Rscript shiny_app/prep/xxx.R`。

#### 2.0 `shiny_app/R/constants.R`（新規）

prep と app の双方が `source()` する定数ファイル。ハードコードをここに集約する。

```r
YEARS      <- c(`1` = 2015L, `2` = 2020L)   # time -> year (FY)
YEAR_RESID <- 2020L                          # 残差・推定値は t=2 のみ
CRITERIA   <- c(pAIC = "pAIC", pBIC = "pBIC")
CRITERION_LABEL <- c(pAIC = "pAIC (penalized AIC)", pBIC = "pBIC (penalized BIC)")
OUTCOMES <- tibble::tribble(
  ~outcome, ~key,             ~label,                                       ~palette,
  "y1",     "aging65",        "Population Aged 65+ (%, standardized)",      "plasma",
  "y2",     "secondary_ind",  "Secondary-Industry Workers (%, standardized)","mako"
)
COVARIATES <- tibble::tribble(
  ~col,            ~key,                 ~label,
  "x_common1",     "ln_taxable_income",  "Log Taxable Income per Taxpayer",
  "x_common2",     "ln_foreign_pop",     "Log(1 + Foreign Population per 100k)",
  "x_specific1_1", "net_inmigration",    "Net In-migration Rate, Japanese migrants (%)",
  "x_specific2_1", "commuter_inflow",    "Commuters from Other Municipalities (%)"
)
MODEL_TYPE / MODEL_FAMILY / MODEL_LABELS / FAMILY_LABEL / FAMILY_COLOR  # 前回 helpers.R から移動
```

#### 2.1 `build_app_data.R`

前回版からの変更点:
- 入力: `output/real/fitted_residuals_<criterion>/`、mapping は `region_id_mapping198.csv`
- `VAR_SPEC` を `constants.R` の `OUTCOMES` / `COVARIATES` から生成
- **`app_data` に列 `criterion`（pAIC / pBIC。応答・説明変数行は NA）を追加**
- 残差行 = `y*_resid`（イノベーション）のみ。`y*_resid_reduced` は読まない。推定値行 = `y*_pred`
- 行の内訳: 応答 2 × 2 年、説明 4 × 2 年、残差 11 モデル × 2 基準 × 2 応答、推定値 同数 → 各 198 行。合計 198 × (4 + 8 + 44 + 44) = **19,800 行**
- `variable` キーは `<outcome_key>_<model>`（例 `aging65_1111`）。criterion は別列
- サニティチェック: 行数 19,800；pAIC と pBIC の残差が 1111/0111/従来 5 モデルで一致し 1011/1101/0101/1001 で不一致；応答の min/max が定義書と一致

#### 2.2 `build_spatial_diag.R`

- ループを `model × outcome × criterion`（44 組）に拡張し、`moran`・`lisa` に `criterion` 列を追加
- 残差はイノベーション `y*_resid`
- `set.seed(20260614)`・`nsim = 999` は前回踏襲
- 期待行数: moran 44、lisa 44 × 198 = 8,712
- **数値ゲート**（pAIC・`moran.test` の I を `reference/heatmaps_Data_I_v2_console.txt` と照合、小数 4 桁）:
  e1(1111) = +0.0169、e2(1111) = +0.0223、e1(dddd) = +0.0073、e2(dddd) = +0.0110

#### 2.3 `build_matrices.R`

- `comparison_table_pAIC.csv` と `_pBIC.csv` の両方を読み、`coef` に `criterion` 列を追加
- 有意性: `signif_psi` が空なら `signif_hessian` を採用し、`signif_source`（"psi" / "hessian"）列を残す
- `fit` 表: `model_selection_summary.csv` に `model_type` を `constants.R` で補い、`pAIC, pBIC, gamma_pAIC, gamma_pBIC, d_eff_pAIC, d_eff_pBIC, delta_pAIC, delta_pBIC` を保持
- `gamma_optimal_summary.csv` を `gamma` として格納 → `coef_matrices.rds = list(coef, fit, gamma)`

#### 2.4 ゲート

3 本すべてが警告なしで完走し、`shiny_app/data/` の 3 つの `.rds` が更新される。各スクリプトのサニティチェック出力を §6 の作業ログに貼る。コミット。

### Phase 3. アプリの改修（1 日）

`app.R` と `helpers.R` は前回版をコピーして差分改修する（5 タブ構成は維持）。

| 箇所 | 変更 |
|---|---|
| 定数 | `helpers.R` の `MODEL_*`、`CATEGORY_LABELS` 等は `constants.R` へ移動し `source` |
| ベースマップ | `add_basemap(map, layerId = NULL, pane = NULL)` ヘルパー。キーがあれば鍵付き URL、無ければ `providers$CartoDB.Positron`。Swipe の左右タイル（pane 指定）も同ヘルパー経由。`stop()` は撤去 |
| 年 | `2016/2021` → `YEARS`（2015/2020）。「2021-only」の文言は `YEAR_RESID` から生成 |
| Tab 1 Explore | `radioButtons("exp_criterion")`（pAIC/pBIC）を追加。`vec_for()` に criterion 引数。残差差分は同一 criterion 内で `e_A − e_B`。**observeEvent のトリガーリストに `input$exp_criterion` を追加**（§5 参照） |
| Tab 2 Panel change | Δ = 2020 − 2015。ラベルは `constants.R` から |
| Tab 3 Swipe | 左右それぞれに model + criterion。既定: 左 0ddd/pAIC、右 0111/pAIC。バッジに criterion と Moran's I（該当 criterion）を表示 |
| Tab 4 Diagnostics | criterion セレクタ追加。Moran 棒グラフは選択中 criterion。任意: pAIC/pBIC を dodge で並べる表示 |
| Tab 5 Model structure | criterion セレクタ追加。fit 表は選択 criterion の pAIC/pBIC でソートし γ, d_eff を表示。subtitle に「Significance: Ψ-based SE; Hessian-based SE where Ψ is unavailable (R, Λ)」 |
| popup | ハードコード（income, pop_growth …）を廃し `OUTCOMES` / `COVARIATES` から生成。残差・推定値表は選択中 criterion の 11 モデル分 |
| 色スケール | 応答 + 推定値は outcome ごと共通範囲（pAIC/pBIC 両方の pred を含める）。残差は全モデル・全基準で対称 ±max。差分は都度 |
| タイトル | `"Kansai MGNST viewer — penalized (pAIC/pBIC)"` |

**ゲート**: launch.json で起動し、5 タブすべてを操作してスクリーンショットを撮る。特に:
- Explore で criterion を pAIC → pBIC に切替えたとき、1011 の残差地図は変わり、1111 は変わらない
- Swipe で左 = 1101 (pAIC)、右 = 1101 (pBIC) が表示でき、左右で色が違う
- キーを外した状態（`Sys.unsetenv`）でもフォールバックで起動する
- R コンソールに warning / error が出ない
- コミット

### Phase 4. 数値照合（半日）

`reference/heatmaps_Data_I_v2_console.txt`（pAIC、モデル 1111 と dddd）と定義書に対し:

| 項目 | 期待値 | 確認方法 |
|---|---|---|
| y1 obs+fitted 共通範囲（1111, dddd のみ） | [-1.889, 3.813] | app_data から同条件で range |
| y2 同上 | [-2.146, 2.995] | 同上 |
| 残差 ±max（1111, dddd のみ） | ±0.930 | 同上（全 11 モデルではこれ以上） |
| Moran's I（§2.2 の 4 値） | 一致 | spatial_diag.rds |
| 記述統計 y1 (2020) mean 0.1752, sd 1.0332 | 一致 | app_data |
| 定義書 §6 の min/max（396 観測） | 一致 | app_data |

### Phase 5. ドキュメント（2 時間）

- `shiny_app/README.md`: 前回 README を土台に、罰則あり版の差分（criterion 軸復活、イノベーション残差、変数・年の変更、real_v2 出典、キーの設定方法とフォールバック）
- `shiny_app/DESIGN_penalized.md`: 本計画 §0〜§3 の確定内容を仕様として整理
- ルート `README.md`: フォルダ構成、再現手順（prep 3 本 → runApp）、**データ出典**（定義書 §8 の引用例をそのまま、N03 行政区域の出典、CARTO/OSM の attribution）、ライセンス（コード MIT、加工データは e-Stat 利用規約に従う）
- `PLAN.md` §6 に作業ログ。コミット

### Phase 6. 環境固定・配布（任意）

- `shiny_app/` で `renv::snapshot()`（`prep/` は `.renvignore` で除外）。ggplot2 4.0.3 など前回 lock より新しい版が入っているため差分は出る
- GitHub へ push（private）
- Posit Connect Cloud へは **新規デプロイ**（前回の `.posit/` は流用しない）。`CARTO_BASEMAP_KEY` を secret に設定。生成された `.posit/publish/*.toml` は前回同様コミットし、`deployments/` は除外

---

## 4. Opus 5.5 への指示のしかた

1. 最初に本 `PLAN.md` と前回の `20260613_noPenalty/shiny_app/`（app.R, helpers.R, prep/, README.md）を読ませる
2. 「§0.1 のソースフォルダは読み取り専用。書き込みは `C:\Users\s.ota\Desktop\残差地図_罰則有\` のみ」を明示する
3. 「`CARTO_BASEMAP_KEY` の値は扱わない。ファイルにも会話にも書かない」を明示する
4. Phase 単位で進めさせ、各ゲートの結果（数値・スクリーンショット）を報告させてから次 Phase を許可する。Phase 2 の Moran's I 照合と Phase 3 の切替動作確認は省略させない
5. R はフルパス `"C:\Program Files\R\R-4.6.1\bin\Rscript.exe"`、コードはファイルに書いて実行（`-e` に長いコードを渡さない）
6. コミットは Phase ごと。push はユーザーが指示したときだけ

---

## 5. 落とし穴リスト（実装者向け）

- `region_id_mapping198_en.csv`（Kansai_Map 側）は列名が違う。使わない
- `output/real/`（旧世代）と `real_v2/` を混ぜない（`pred` の定義が違う）
- `y*_resid_reduced` 列は読まない。うっかり `resid` と取り違えると空間モデルの残差 Moran's I が変わる
- 残差 CSV の `time` は 2 のみ（198 行）。`filter(time == 2)` は残すが、行数チェックで 198 を要求する
- `Sigma` は `[1,1],[1,2],[2,2]` のみ格納。前回どおり `[2,1]` をミラーする
- `A` 行は時間ラグ無しモデル（1101, 0101, 1001, 000d）で NA → `present = FALSE`
- `R[1,2]`, `R[2,1]`, `Lambda[1,2]`, `Lambda[2,1]`, `Sigma[1,2]` が従来モデルで厳密に 0 → 構造的 0（灰色）。前回ロジックを維持
- 前回 `app.R` の Explore タブは `observeEvent` のトリガーリストに読む入力を全部列挙している（依存追跡の不具合回避）。criterion を追加したら **トリガーリストにも必ず追加**
- leaflet.extras2 の side-by-side は「タイルレイヤ」の layerId で左右を判定する。前回のペイン構成（`left` / `right` pane にタイルとポリゴンを同居）を崩さない。`add_basemap()` ヘルパーでもこの構成を保つ
- `x_specific2_1`（通勤者比率）は最大 231 と分布が偏る。線形 viridis で問題ないが、必要なら分位分類の追加を検討（今回は必須にしない）
- Google Drive（G:）上で sf を使う R を直接走らせない。C: の本フォルダで作業する
- 秘密情報: `.Renviron` と `.posit/publish/deployments/` は `.gitignore` に入れ、`git status` で追跡されていないことを最初のコミット前に確認する

---

## 6. 作業ログ（Opus が追記）

（Phase ごとの実行コマンド、ゲート結果、判断を追記する）

### Phase 0（2026-09-23, Opus 5.5）

- リポジトリ名: **`mstr-RSS-2026-kansai-penalized-map`**（ユーザー指定）
- `.claude/launch.json` 作成（R 4.6.1, port 7897, name `shiny-penalized`）
- `.gitignore`: 前回から流用、`shiny_app_alpha/` 行を削除、`.Renviron` を追加
- `git init -b main`、`git remote add origin https://github.com/s-ohta-s/mstr-RSS-2026-kansai-penalized-map.git`（未 push。GitHub 側の空 private リポジトリはユーザーが作成）
- ユーザーが `Documents\.Renviron` を作成（BOM なし、末尾改行あり。値は扱っていない）
  - 初回は値がプレースホルダ `<キー>` のままで、`!= ''` だけのゲートを素通りした → 形の検査（`<>`・非 ASCII・前後空白・引用符なし）を追加し、ユーザー修正後に全項目 OK を確認
  - キーが CARTO で実際に有効かは Phase 3 のタイル表示で確認する
- **ゲート通過**: R 4.6.1、shiny 1.14.0 / leaflet 2.2.3 / leaflet.extras2 1.3.2 / sf 1.1.1 / dplyr 1.2.1 / ggplot2 4.0.3 / spdep 1.4.2 / viridisLite 0.4.3 / RColorBrewer 1.1.3 / renv 1.2.4 / htmltools 0.5.9、`CARTO_BASEMAP_KEY set: TRUE`
- 判断: Git Bash から Rscript を呼ぶと `HOME=C:\Users\s.ota` が引き継がれ `Documents\.Renviron` が読まれない。**R は PowerShell から `Remove-Item Env:HOME` 後に実行する**（launch.json 起動と同条件）
- コミットは PLAN どおり Phase 1 で初回

### Phase 1（2026-09-23, Opus 5.5）

- PowerShell `Copy-Item` で 38 ファイルをコピー（§2 の一覧どおり。`real_v2/` の xlsx・`gamma_search_*`・`se_comparison_*`・`identifiability_*` はコピーしない）。**SHA256 はコピー元と全件一致（不一致 0）**
- ゲート（`Rscript` スクリプト、作業ディレクトリ = 本フォルダ）: **PASS**
  - fitted_residuals: pAIC 11 本 / pBIC 11 本、各 198 行・10 列・`time == 2`・region 1..198・NA なし
  - **列順は §0.2 の記述と異なる**: 実際は `region,time,y1_obs,y1_pred,y1_resid,y1_resid_reduced,y2_obs,y2_pred,y2_resid,y2_resid_reduced`（列の集合は同じ）。prep は列名で参照するので影響なし
  - mapping: 列 `region,code,name`、198 行
  - transformed: 396 行（time 1/2 × 198）、min/max は定義書と全 6 列一致（y1 [-1.8887, 3.5265]、y2 [-2.1456, 2.9951]、x_common1 [7.7212, 8.7828]、x_common2 [0, 8.6415]、x_specific1_1 [-5.79, 1.41]、x_specific2_1 [7.8, 231.1]）
  - muni_sf: is_target TRUE 198 / FALSE 7、対象コードと `substr(mapping$code,1,5)` が 198 件全一致、EPSG 4326。pref_borders 13 行
  - W: 198×198、行和 = 1
- コミット前確認: 改行はすべて LF（`core.autocrlf=true` でも内容は変わらない）。`.claude/` は ignored。秘密情報らしき文字列の検索でヒットしたのは renv.lock のパッケージ説明文と、変数名に触れた `.gitignore`/PLAN.md のみ
