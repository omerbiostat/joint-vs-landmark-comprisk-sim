# =====================================================================
#  JM vs LM  —  PUBLICATION FIGURES  (from the CSV run, w=6 excluded)
#  Paketler:  install.packages(c("ggplot2","dplyr","tidyr","patchwork"))
#  Run:  source("figures.R")
#  Output: Figure1, Figure2, FigureS1, FigureS2  (300 dpi TIFF + PNG)
# =====================================================================
library(ggplot2); library(dplyr); library(tidyr); library(patchwork)

data_dir <- "path/to/results-folder"
files <- c(
  "s2.5 / w2" = "results_s2.5_w2.csv",
  "s2.5 / w4" = "results_s2.5_w4.csv",
  "s5.5 / w2" = "results_s5.5_w2.csv",
  "s5.5 / w4" = "results_s5.5_w4.csv"
)

## --- read all replications into one long table ---
raw <- bind_rows(lapply(names(files), function(sc) {
  d <- read.csv(file.path(data_dir, files[[sc]]))
  d$scenario <- sc; d
}))
raw$scenario <- factor(raw$scenario, levels = names(files))

pal <- c("JM" = "#D55E00", "LM" = "#0072B2")          # colour-blind safe
theme_set(theme_bw(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  strip.background = element_rect(fill = "grey92", color = NA),
                  legend.position  = "right"))

## =====================================================================
## FIGURE 1 — Performance (JM vs LM) + Monte Carlo standard error (MCSE)
## =====================================================================

grid <- expand.grid(metric = c("bias","rmse","brier","auc"),
                    outcome = c("trans","death"),
                    method  = c("jm","lm"), stringsAsFactors = FALSE)
long1 <- bind_rows(lapply(seq_len(nrow(grid)), function(k) {
  col <- paste(grid$metric[k], grid$outcome[k], grid$method[k], sep = "_")
  raw %>% group_by(scenario) %>%
    summarise(mean = mean(.data[[col]], na.rm = TRUE),
              mcse = sd(.data[[col]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[col]]))),
              .groups = "drop") %>%
    mutate(metric = grid$metric[k], outcome = grid$outcome[k],
           method = toupper(grid$method[k]))
}))
long1$metric  <- factor(long1$metric,  levels = c("bias","rmse","brier","auc"),
                        labels = c("Bias","RMSE","Brier","AUC"))
long1$outcome <- factor(long1$outcome, levels = c("trans","death"),
                        labels = c("Transplantation","Death"))

fig1 <- ggplot(long1, aes(scenario, mean, color = method, group = method)) +
  geom_hline(data = subset(long1, metric == "Bias"),
             aes(yintercept = 0), linetype = "dotted", color = "grey55") +
  
  # 0.50 reference line for AUC (no label)
  geom_hline(data = subset(long1, metric == "AUC"),
             aes(yintercept = 0.50), linetype = "dashed", color = "grey55") +
  
  geom_errorbar(aes(ymin = mean - mcse, ymax = mean + mcse),
                width = .18, position = position_dodge(.45)) +
  geom_point(size = 2.2, position = position_dodge(.45)) +
  facet_grid(metric ~ outcome, scales = "free_y", switch = "y") +
  scale_color_manual(values = pal, name = NULL) +
  labs(x = NULL, y = NULL) + # ,caption = "Points: mean over 200 replications. Error bars: Monte Carlo standard error."
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        strip.placement = "outside")
ggsave("Figure1_performance.tiff", fig1, width = 7.0, height = 8.2, dpi = 300, compression = "lzw")
ggsave("Figure1_performance.png",  fig1, width = 7.0, height = 8.2, dpi = 300)

## =====================================================================
## FIGURE 2 — MECHANISM (A: redistribution · B: low-EPV instability)
## =====================================================================
# 2A — bias_death ~ bias_trans, y = -x line, r per panel
sc2a <- bind_rows(
  raw %>% transmute(scenario, method = "JM", bt = bias_trans_jm, bd = bias_death_jm),
  raw %>% transmute(scenario, method = "LM", bt = bias_trans_lm, bd = bias_death_lm)
)
rlab <- sc2a %>% group_by(scenario, method) %>%
  summarise(r = cor(bt, bd, use = "complete.obs"), .groups = "drop") %>%
  mutate(t = sprintf("%s: r = %.2f", method, r)) %>%
  group_by(scenario) %>% summarise(lab = paste(t, collapse = "\n"), .groups = "drop")

fig2a <- ggplot(sc2a, aes(bt, bd, color = method)) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, color = "grey88") +
  geom_hline(yintercept = 0, color = "grey88") +
  geom_point(alpha = .55, size = 1.2) +
  geom_text(data = rlab, aes(x = Inf, y = Inf, label = lab),
            hjust = 1.05, vjust = 1.12, size = 3.0, color = "black", inherit.aes = FALSE) +
  facet_wrap(~ scenario, nrow = 1) +
  scale_color_manual(values = pal, name = NULL) +
  labs(x = "Prediction error — transplantation (bias)",
       y = "Prediction error — death (bias)",
       subtitle = "Redistribution of probability mass between competing causes (dashed line: y = -x)")  + coord_fixed(ratio = 1)

# 2B — transplant bias density: JM heavy right tail vs LM tight at 0
sc2b <- bind_rows(
  raw %>% transmute(scenario, method = "JM", bt = bias_trans_jm),
  raw %>% transmute(scenario, method = "LM", bt = bias_trans_lm)
)
fig2b <- ggplot(sc2b, aes(bt, fill = method, color = method)) +
  geom_density(alpha = .35) +
  geom_vline(xintercept = 0, color = "grey60", linetype = "dotted") +
  facet_wrap(~ scenario, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = pal, name = NULL) +
  scale_color_manual(values = pal, name = NULL) +
  labs(x = "Per-replication bias — transplantation", y = "Density",
       subtitle = "Joint Model bias arises from a minority of unstable, low-EPV replications (heavy right tail)") # + guides(fill = "none", color = "none")

fig2 <- fig2a / fig2b + plot_layout(guides = "collect", heights = c(1.40, 1)) + plot_annotation(tag_levels = 'A')  # & theme(legend.position = "right") 
fig2
ggsave("Figure2_mechanism.tiff", fig2, width = 9.5, height = 7.0, dpi = 300, compression = "lzw")
ggsave("Figure2_mechanism.png",  fig2, width = 9.5, height = 7.0, dpi = 300)

## =====================================================================
## FIGURE S1 — Coverage (JM vs LM), nominal 0.95
## =====================================================================
covsum <- bind_rows(
  raw %>% group_by(scenario) %>%
    summarise(method = "JM",
              Transplantation = mean(cov_trans_jm, na.rm = TRUE),
              Death           = mean(cov_death_jm, na.rm = TRUE), .groups = "drop"),
  raw %>% group_by(scenario) %>%
    summarise(method = "LM",
              Transplantation = mean(cov_trans_lm, na.rm = TRUE),
              Death           = mean(cov_death_lm, na.rm = TRUE), .groups = "drop")
) %>% pivot_longer(c(Transplantation, Death), names_to = "outcome", values_to = "coverage")
covsum$outcome <- factor(covsum$outcome, levels = c("Transplantation","Death"))

figS1 <- ggplot(covsum, aes(scenario, coverage, fill = method)) +
  geom_col(position = position_dodge(.7), width = .6) +
  geom_hline(yintercept = .95, linetype = "dashed", color = "grey40") +
  facet_wrap(~ outcome) +
  scale_fill_manual(values = pal, name = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Coverage of 95% intervals") + # ,caption = "Dashed line: nominal 0.95. LM from bootstrap on a subset of replications."
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave("FigureS1_coverage.tiff", figS1, width = 7.5, height = 4.0, dpi = 300, compression = "lzw")
ggsave("FigureS1_coverage.png",  figS1, width = 7.5, height = 4.0, dpi = 300)

## =====================================================================
## FIGURE S2 — Per-replication bias distributions (all scenarios)
## =====================================================================
scS2 <- bind_rows(
  raw %>% transmute(scenario, method = "JM",
                    Transplantation = bias_trans_jm, Death = bias_death_jm),
  raw %>% transmute(scenario, method = "LM",
                    Transplantation = bias_trans_lm, Death = bias_death_lm)
) %>% pivot_longer(c(Transplantation, Death), names_to = "outcome", values_to = "bias")
scS2$outcome <- factor(scS2$outcome, levels = c("Transplantation","Death"))

figS2 <- ggplot(scS2, aes(scenario, bias, fill = method)) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dotted") +
  geom_violin(position = position_dodge(.8), width = .9, alpha = .4, color = NA) +
  geom_boxplot(position = position_dodge(.8), width = .18, outlier.size = .4, alpha = .9, outlier.alpha = 0.4) +
  facet_wrap(~ outcome) + # , scales = "free_y"
  coord_cartesian(ylim = c(-0.30, 0.30)) +
  scale_fill_manual(values = pal, name = NULL) +
  labs(x = NULL, y = "Per-replication bias") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave("FigureS2_bias_distributions.tiff", figS2, width = 8.5, height = 5.5, dpi = 300, compression = "lzw")
ggsave("FigureS2_bias_distributions.png",  figS2, width = 8.5, height = 5.5, dpi = 300)

cat("\n[OK] Figures written: Figure1, Figure2, FigureS1, FigureS2 (TIFF + PNG).\n")
cat("PLOS requires TIFF/EPS; TIFF files are 300 dpi. Adjust sizes via ggsave width/height if needed.\n")
