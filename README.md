# 都市化勾配 × カエル類メタ解析（再現リポジトリ）

このリポジトリは、都市化勾配に対するカエル類の **Richness / Abundance** のメタ解析データ（CSV）と、
解析・図を再現する **Rスクリプト**を公開します。誌面掲載図はフォントや余白を微調整する場合がありますが、
**推定値・信頼区間は同一**です。

This repo provides the meta-analysis dataset (CSV) and R scripts to reproduce analyses and base figures.
Published figures may have minor layout tweaks (fonts/margins), but **estimates and CIs are identical**.

---

## 実行方法
```r
# 依存パッケージ（初回のみインストール例）
# install.packages(c("readr","dplyr","metafor","ggplot2","ggrepel"))

# 実行
source("R/run_meta.R")
```
- 出力先: `data/ForAnalysis/figures/`（PDF）
- バブル図の重み: **面積 ∝ 1/SE²**（逆分散重み）、**半径 ∝ 1/SE**
- データ読み込みは相対パス（または環境変数 `META_DATA_PATH`）を使用

---

## ファイル
- `data/ForAnalysis/meta_effects_master_20250705_wPop_cleaned.csv`（入力データ）
- `R/run_meta.R`（再現スクリプト：フォレスト／バブル／勾配別の図を生成）

---

## 再現性メモ
- 依存関係を固定したい場合（任意）：
```r
install.packages("renv")
renv::init()
renv::snapshot()   # → 生成された renv.lock をコミット
# 別環境では renv::restore() で再現
```

---

## ライセンス
- **Code:** MIT  
- **Data:** CC BY 4.0

## 引用例（任意）
> Data and code to reproduce the meta-analyses and base figures are available at: ＜本リポのURLまたはDOI＞.
