# 設計ドキュメント（罰則あり版）: MGNST 結果可視化 Shiny アプリ — penalized (pAIC/pBIC)

前回（罰則なし版）の `DESIGN.md` / `DESIGN_beta.md`（`20260613_noPenalty/shiny_app/`）からの**差分仕様**。
書かれていない点は前回の設計を踏襲する。作業の経緯と検証結果は、リポジトリ直下の `PLAN.md` §6（作業ログ）にある。

---

## 1. 目的と基本方針（前回から変更なし）

- テーマ: 新しい多変量空間手法（MGNS 系）と、従来の個別回帰モデルを**地図と診断で視覚的に比較**する
- **アプリは読むだけ**。推定・Moran's I・LISA・行列の整形はすべて `prep/` で事前計算して `.rds` に落とす。アプリ実行時に `spdep` は使わない
- 5 タブ構成（Explore / Panel change / Compare (swipe) / Spatial diagnostics / Model structure）を維持
- UI は英語。タイトル `Kansai MGNST viewer — penalized (pAIC/pBIC)`

## 2. 入力（すべて本リポジトリ内のコピー）

| 種類 | ファイル | 元の場所 |
|---|---|---|
| 変数 | `data/real/transformed_data198_for_R.csv`（396 行） | `20260705_Kansai_Estimation/data/real/` |
| 対応表 | `data/real/region_id_mapping198.csv`（列 `region, code, name`） | 同上 |
| 空間重み | `data/real/W_row_standardized_198_for_R.csv` | 同上 |
| 定義書 | `data/real/DATA_SOURCE_AND_PROCESSING.md`, `data_source_and_processing_en.md` | `20260702_RealDataKansai/Data_I/` |
| 推定結果 | `output/real/`（comparison_table_pAIC/pBIC, model_selection_summary, gamma_optimal_summary, fitted_residuals_pAIC/pBIC） | `20260705_Kansai_Estimation/real_v2/` |
| 照合用 | `reference/create_heatmaps_en.r`, `heatmaps_Data_I_v2_console.txt` | `20260709_Kansai_Map/` |
| 地理 | `shiny_app/data/muni_sf.rds`, `pref_borders.rds` | 前回アプリから流用（推定に依存しない） |

- 推定結果は **`real_v2/`** のみを使う。旧世代 `output/real/`（8 列）とは `pred` の定義が違うので混ぜない
- 残差 CSV の列（実際の順）: `region, time, y1_obs, y1_pred, y1_resid, y1_resid_reduced, y2_obs, y2_pred, y2_resid, y2_resid_reduced`。`time == 2` の 198 行のみ

## 3. 変数・年・基準（`R/constants.R` に一元化）

prep と app の両方が `source` する。ハードコードはここ以外に置かない。

| 定数 | 内容 |
|---|---|
| `YEARS` | `time` 1 → 2015、2 → 2020（FY） |
| `YEAR_BASE` / `YEAR_RESID` | 2015 / 2020（推定値・残差は 2020 のみ） |
| `CRITERIA`, `CRITERION_LABEL` | `pAIC`, `pBIC` |
| `OUTCOMES` | y1 = `aging65` "Population Aged 65+ (%, standardized)"（plasma）、y2 = `secondary_ind` "Secondary-Industry Workers (%, standardized)"（mako）。短縮名 "Aged 65+" / "Secondary industry" |
| `COVARIATES` | `x_common1` = `ln_taxable_income`、`x_common2` = `ln_foreign_pop`、`x_specific1_1` = `net_inmigration`、`x_specific2_1` = `commuter_inflow` |
| `CATEGORY_LABELS`, `RESID_CATS` | Explore のカテゴリ。`residual`, `fitted`, `residual_diff` は 2020 のみ・criterion 依存 |
| `MODEL_TYPE`, `MODEL_IDS`, `MODEL_LABELS`, `MODEL_FAMILY`, `FAMILY_LABEL`, `FAMILY_COLOR` | 11 モデルの表示名と比較軸（new / conventional） |

## 4. 残差と推定値の定義

- **残差 = イノベーション残差 `y*_resid`**（尤度・Σ 推定の残差）。地図、Moran's I、LISA はすべてこれを使う
- `y*_resid_reduced`（= obs − pred、誘導形残差）は**読まない**。residual-type の切替は作らない
- **推定値 = `y*_pred`**（誘導形トレンド予測 ŷ。OLS/VARX/MSEM は Xβ̂、MSAR/MGNS は (I − R̂⊗W)⁻¹Xβ̂）

## 5. criterion（pAIC / pBIC）軸

- 6 つの罰則モデルは基準ごとに最適 γ を持つ。本推定では **1011, 1101, 0101, 1001 は pAIC と pBIC で残差・推定値・行列が異なり**、1111, 0111（γ が同じ）と従来 5 モデル（γ = 0）は同一
- UI:
  - Explore: criterion ラジオ（残差・推定値・残差差分のときだけ表示）。残差差分は同一 criterion 内で `e_A − e_B`
  - Swipe: 左右それぞれに model と criterion（同じモデルの pAIC vs pBIC を比較できる）。既定は左 0ddd/pAIC、右 0111/pAIC
  - Diagnostics / Model structure: criterion ラジオ
- Explore の `observeEvent` はトリガーリストに読む入力をすべて列挙する（依存追跡の不具合回避）。**`input$exp_criterion` を含めること**

## 6. データ構造（`.rds`）

### 6.1 `app_data.rds`（19,800 行）

列: `muni_code, year, category, variable, outcome, model, criterion, value, label`

- 内訳: 応答 2 × 2 年、説明変数 4 × 2 年、残差 11 モデル × 2 基準 × 2 応答、推定値も同数（各 198 行）
- `variable` キー: 応答・説明変数は `OUTCOMES$key` / `COVARIATES$key`、残差・推定値は `<outcome key>_<model>`（例 `aging65_1111`）。**criterion は別列**で、キーと label は pAIC/pBIC で共通（選択肢が重複しない）
- label 例: `Residual: Aged 65+ — MGNS (1111)`

### 6.2 `spatial_diag.rds` = `list(moran, lisa)`

- `moran`（44 行 = 11 モデル × 2 応答 × 2 基準）: `morans_i, expectation, variance, z, p_analytic`（`moran.test`, randomisation）、`p_value`（`moran.mc`, 999 回, `set.seed(20260614)`）
- `lisa`（8,712 行）: `local_i, quadrant`（HH/LL/HL/LH）、`p_value`（`localmoran_perm`）。NS は app で p < α から導く

### 6.3 `coef_matrices.rds` = `list(coef, fit, gamma)`

- `coef`（352 行 = (15 セル + Σ[2,1] のミラー) × 11 × 2）: `estimate, se, signif, signif_source, present`
  - `present`: NA → モデルに無い、非対角の厳密な 0 → 構造的 0（灰色）
  - **有意性**: Ψ-SE があれば `signif_psi`、無く Hessian-SE があれば `signif_hessian`、どちらも無ければ記号なし。`signif_source` = `"psi"` / `"hessian"` / NA。判定は**SE の有無**で行う（`signif_psi` が空欄でも Ψ-SE があれば「有意でない」の意味なので、Hessian に切り替えない）
  - 実際の割当: A = Ψ、R と Λ（新手法）= Hessian、Σ と従来モデルの R/Λ = なし
- `fit`（11 行）: `model_selection_summary.csv` + `MODEL_TYPE`/`MODEL_FAMILY`（元 CSV に `model_type` 列が無い）。`d_eff_*, loglik_*, pAIC, pBIC, gamma_*, delta_*` を保持
- `gamma`（6 行）: `gamma_optimal_summary.csv`
- `model_id` は `colClasses = "character"` で読む（数字だけの ID が数値化され `0111` → `111` になるのを防ぐ）

## 7. 色スケール

- パレットは `reference/create_heatmaps_en.r` を踏襲: y1 plasma、y2 mako、説明変数 viridis、残差・差分 RdBu（`direction = 1` 相当、低 = 赤・高 = 青）、Δ も RdBu 0 中心
- 応答 + 推定値: outcome ごとに 1 つの範囲（観測値 2 年 + 全 11 モデル × 2 基準の推定値）
- 残差: 全モデル・全応答・全基準で対称 ±max
- 残差差分・Δ・Swipe: 表示中のデータから都度
- 凡例は小さい値が下（反転パレット + ラベル降順）
- **既知の性質**: 共通範囲は当てはまりの悪い 000d で決まる（y1 [-3.03, 3.89]、残差 ±2.69）。良いモデルの残差（±0.93 以内）は Explore で薄い色になる。今回は変更しない（代案は `PLAN.md` Phase 4 のログ）

## 8. ベースマップ

- `add_basemap(map, layerId = NULL, pane = NULL)`: `CARTO_BASEMAP_KEY` があれば鍵付き CARTO Positron URL、無ければ `providers$CartoDB.Positron`
- Swipe（leaflet.extras2 side-by-side）は「タイル層」を layerId で探し、pane で切り抜く。左右それぞれ `left` / `right` pane にタイル 1 枚とポリゴンを同居させる構成を維持（`addProviderTiles` も tile として登録されるので、フォールバック時も同じ構成で動く）
- 起動時に使用中のベースマップを `message()` で出す（キーの値は出さない）
- **キーなしの場合、CARTO は「API KEY REQUIRED」の透かし入りタイルを返す**。アプリは動く。代替プロバイダには変えない（2026-09-23 ユーザー判断）
- キーはリポジトリに入れない（`.gitignore` に `.Renviron`）。ローカルは `Documents\.Renviron`、Posit Connect Cloud は secret

## 9. ホバー / クリック

- ホバー: 市町村名（英・日）+ 表示中の値。ラベルは **名前なしのリスト**で渡す（`mapply(..., USE.NAMES = FALSE)`）。名前付きだと JSON オブジェクトになり、leaflet が全ツールチップに 198 件分を渡して表示に失敗する（前回アプリにあった不具合）
- クリック: `popup_all()` で 198 件を一括生成。選択年の応答・説明変数と、選択 criterion の 11 モデルの残差・推定値（y1 | y2）

## 10. 検証（PLAN.md §6 に詳細）

- prep 3 本: 警告 0、サニティチェック（行数、pAIC/pBIC の一致・不一致パターン、定義書の範囲、muni_sf との整合）
- Moran's I（pAIC, 1111/dddd）が `reference/heatmaps_Data_I_v2_console.txt` と一致
- アプリ: criterion 切替で 1011 の残差地図が変わり 1111 は変わらない、Swipe で 1101 pAIC vs pBIC、キーなし起動
- 数値照合 38 項目（統計表 20 行、共通範囲、Moran's I の I/SD/z/p、残差差分、定義書 §6）がすべて一致
