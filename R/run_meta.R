rm(list = ls())  # すべてのオブジェクトを削除

# 必要パッケージ
library(readr); library(dplyr); library(metafor)
library(ggplot2); library(ggrepel)

# ==== 入力（相対パス or 環境変数）====
data_path <- Sys.getenv(
  "META_DATA_PATH",
  unset = file.path("data","ForAnalysis","meta_effects_master_20250705_wPop_cleaned.csv")
)
meta <- read_csv(data_path, show_col_types = FALSE)

# 共通ユーティリティ
prep_subset <- function(df) {
  df %>%
    distinct(StudyID, .keep_all = TRUE) %>%      # ID重複防止（先頭行を採用）
    filter(!is.na(z), !is.na(SE_z), SE_z > 0)    # rma()の安全網
}
date_tag <- format(Sys.Date(), "%Y%m%d")
out_dir  <- "data/ForAnalysis/figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Richness (Full)：フォレスト ─────────────────────────────
dat_full <- meta %>%
  filter(Outcome == "Richness", Gradient == "Full") %>%
  prep_subset()

res_full <- rma(yi = z, sei = SE_z, data = dat_full, method = "REML")

forest(res_full,
       transf = transf.ztor,                   # z→r を明示
       xlab   = "Correlation (r)",
       slab   = gsub("_", "", dat_full$AuthorYear, fixed = TRUE),
       main   = "Richness (Full)")

pdf(file.path(out_dir, paste0("forest_richness_full_", date_tag, ".pdf")), width = 7, height = 5)
forest(res_full,
       transf = transf.ztor,
       xlab   = "Correlation (r)",
       slab   = gsub("_", "", dat_full$AuthorYear, fixed = TRUE),
       main   = "Richness (Full)")
dev.off()

# ── Richness (Full) × LogPop：バブル（面積 ∝ 1/SE^2） ───────
dat_city <- dat_full %>%
  mutate(r = transf.ztor(z),
         weight = 1 / (SE_z^2)) %>%
  filter(is.finite(LogPop))

# メタ回帰（数値はrmaで推定。geom_smoothは視覚トレンド）
res_reg <- rma(yi = z, sei = SE_z, mods = ~ LogPop, data = dat_city, method = "REML")

p <- ggplot(dat_city, aes(x = LogPop, y = r)) +
  geom_point(aes(size = weight), color = "black", alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 1) +
  geom_text_repel(aes(label = gsub("_", "", AuthorYear, fixed = TRUE)),
                  size = 3.5, force = 1, nudge_y = 0.03,
                  segment.color = NA, max.overlaps = 10, box.padding = 0.3) +
  scale_size_area(name = "Inverse-variance weight (1/SE²)", max_size = 10) +
  theme_minimal(base_size = 14) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
        panel.grid.major = element_line(color = "gray85", linewidth = 0.4),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(title = "Richness (Full × LogPop)",
       x = "ln(city population)", y = "Correlation (r)")
print(p)

ggsave(file.path(out_dir, paste0("bubble_richness_logpop_", date_tag, ".pdf")),
       width = 7, height = 4.5, dpi = 300)

# ── Abundance (Full)：フォレスト ────────────────────────────
dat_abundance_full <- meta %>%
  filter(Outcome == "Abundance", Gradient == "Full") %>%
  prep_subset()

res_abundance_full <- rma(yi = z, sei = SE_z, data = dat_abundance_full, method = "REML")

forest(res_abundance_full,
       transf = transf.ztor,
       xlab   = "Correlation (r)",
       slab   = gsub("_", "", dat_abundance_full$AuthorYear, fixed = TRUE),
       main   = "Abundance (Full)")

pdf(file.path(out_dir, paste0("forest_abundance_full_", date_tag, ".pdf")), width = 9, height = 5)
forest(res_abundance_full,
       transf = transf.ztor,
       xlab   = "Correlation (r)",
       slab   = gsub("_", "", dat_abundance_full$AuthorYear, fixed = TRUE),
       main   = "Abundance (Full)")
dev.off()

# ── Richness：勾配別（Urban–Rural系）────────────────────────
# 1) Urban–Rural（Rice除外）
plot_data <- meta %>%
  filter(Outcome == "Richness", Gradient == "Urban_Rural") %>%
  prep_subset()

res_urban_rural_norice <- rma(yi = z, sei = SE_z, data = plot_data, method = "REML")
forest(res_urban_rural_norice,
       transf = transf.ztor,
       xlab = "Correlation (r)",
       slab = gsub("_", "", plot_data$AuthorYear, fixed = TRUE),
       main = "Richness (Urban–Rural Gradient without Paddy Field)")

pdf(file.path(out_dir, paste0("forest_richness_urban_rural_without_rice_", date_tag, ".pdf")), width = 7.5, height = 4.5)
forest(res_urban_rural_norice,
       transf = transf.ztor,
       xlab   = "Correlation (r)",
       slab   = gsub("_", "", plot_data$AuthorYear, fixed = TRUE),
       main   = "Richness (Urban–Rural Gradient without Paddy Field)")
dev.off()

# 2) Urban–Rural（Rice含む：Urban_Rural + Urban_Rice）
dat_urban_rural_withrice <- meta %>%
  filter(Outcome == "Richness", Gradient %in% c("Urban_Rural","Urban_Rice")) %>%
  prep_subset()

res_urban_rural_withrice <- rma(yi = z, sei = SE_z, data = dat_urban_rural_withrice, method = "REML")
forest(res_urban_rural_withrice,
       transf = transf.ztor,
       xlab   = "Correlation (r)",
       slab   = gsub("_", "", dat_urban_rural_withrice$AuthorYear, fixed = TRUE),
       main   = "Richness (Urban–Rural Gradient with Paddy Field)")

pdf(file.path(out_dir, paste0("forest_richness_urban_rural_with_rice_", date_tag, ".pdf")), width = 7.5, height = 4.5)
forest(res_urban_rural_withrice,
       transf = transf.ztor,
       xlab   = "Correlation (r)",
       slab   = gsub("_", "", dat_urban_rural_withrice$AuthorYear, fixed = TRUE),
       main   = "Richness (Urban–Rural Gradient with Paddy Field)")
dev.off()

# 3) Urban–Suburban
dat_us <- meta %>%
  filter(Outcome == "Richness", Gradient == "Urban_Suburban") %>%
  prep_subset()

res_us <- rma(yi = z, sei = SE_z, data = dat_us, method = "REML")
forest(res_us,
       transf = transf.ztor,
       xlab = "Correlation (r)",
       slab = gsub("_", "", dat_us$AuthorYear, fixed = TRUE),
       main = "Richness (Urban–Suburban Gradient)")

pdf(file.path(out_dir, paste0("forest_richness_urban_suburban_", date_tag, ".pdf")), width = 9, height = 5)
forest(res_us,
       transf = transf.ztor,
       xlab = "Correlation (r)",
       slab = gsub("_", "", dat_us$AuthorYear, fixed = TRUE),
       main = "Richness (Urban–Suburban Gradient)")
dev.off()

# 4) Suburban–Rural（Rice除外／含む）
dat_sr_without_rice <- meta %>%
  filter(Outcome == "Richness", Gradient == "Suburban_Rural") %>%
  prep_subset()

res_sr_norice <- rma(yi = z, sei = SE_z, data = dat_sr_without_rice, method = "REML")
forest(res_sr_norice,
       transf = transf.ztor,
       xlab = "Correlation (r)",
       slab = gsub("_", "", dat_sr_without_rice$AuthorYear, fixed = TRUE),
       main = "Richness (Suburban–Rural Gradient without Paddy Field)")

pdf(file.path(out_dir, paste0("forest_richness_suburban_rural_without_rice_", date_tag, ".pdf")), width = 9, height = 5)
forest(res_sr_norice,
       transf = transf.ztor,
       xlab = "Correlation (r)",
       slab = gsub("_", "", dat_sr_without_rice$AuthorYear, fixed = TRUE),
       main = "Richness (Suburban–Rural Gradient without Paddy Field)")
dev.off()

dat_sr_with_rice <- meta %>%
  filter(Outcome == "Richness", Gradient %in% c("Suburban_Rural", "Suburban_Rice")) %>%
  prep_subset()

res_sr_withrice <- rma(yi = z, sei = SE_z, data = dat_sr_with_rice, method = "REML")
forest(res_sr_withrice,
       transf = transf.ztor,
       xlab = "Correlation (r)",
       slab = gsub("_", "", dat_sr_with_rice$AuthorYear, fixed = TRUE),
       main = "Richness (Suburban–Rural Gradient with Paddy Field)")

pdf(file.path(out_dir, paste0("forest_richness_suburban_rural_with_rice_", date_tag, ".pdf")), width = 9, height = 5)
forest(res_sr_withrice,
       transf = transf.ztor,
       xlab = "Correlation (r)",
       slab = gsub("_", "", dat_sr_with_rice$AuthorYear, fixed = TRUE),
       main = "Richness (Suburban–Rural Gradient with Paddy Field)")
dev.off()

# 解析環境ログ
writeLines(capture.output(sessionInfo()), file.path(out_dir, paste0("sessionInfo_meta_", date_tag, ".txt")))
