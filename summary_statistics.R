rm(list = ls()) # .rs.restartR()
Sys.setlocale(category = "LC_ALL", locale = "pt_PT.UTF-8")

# =====================================================================
#  JM vs LM  —  200-REPLICATION SUMMARY STATISTICS
#  Four CSVs (n=1000, w=2,4)  ->  console summary + publication-ready table CSVs
#  Base R only (no extra packages required).
#  Run:  source("summary_statistics.R")
# =====================================================================
setwd("path/to/results-folder")
getwd()
## ---- 1) SETTINGS --------------------------------------------------------
data_dir <- "path/to/results-folder"   # folder containing the CSVs (write the full path if needed)

files <- list(
  "s2.5_w2" = "results_s2.5_w2.csv",
  "s2.5_w4" = "results_s2.5_w4.csv",
  "s5.5_w2" = "results_s5.5_w2.csv",
  "s5.5_w4" = "results_s5.5_w4.csv"
)
scen_s <- c(2.5, 2.5, 5.5, 5.5)
scen_w <- c(2,   4,   2,   4)

# Mechanism — EPV denominator: number of parameters in the JM survival sub-model.
# !!! Verify against the ACTUAL number in your own model. Default ~28.
n_params_surv <- 28

# "Breaking" replication threshold (on transplant bias):
break_thr <- 0.05

## ---- 2) Helpers -------------------------------------------------
num  <- function(x) suppressWarnings(as.numeric(x))
msd  <- function(x) {                      # mean, sd, mcse, n (NA-robust)
  x <- num(x); k <- sum(!is.na(x))
  c(mean = mean(x, na.rm = TRUE),
    sd   = sd(x,   na.rm = TRUE),
    mcse = sd(x,   na.rm = TRUE) / sqrt(k),
    n    = k)
}
fmt  <- function(m, s) sprintf("%.3f (%.3f)", m, s)   # "mean (sd)"

metrics  <- c("bias", "rmse", "brier", "auc")
outcomes <- c("trans", "death")
methods  <- c("jm", "lm")

perf_rows <- mcse_rows <- cov_rows <- mech_rows <- list()

## ---- 3) Compute for each scenario -----------------------------------
for (i in seq_along(files)) {
  sc <- names(files)[i]
  df <- read.csv(file.path(data_dir, files[[sc]]))   # "NA" -> NA automatically
  s_i <- scen_s[i]; w_i <- scen_w[i]

  ## 3a) Performance (mean (SD))  + MCSE
  perf <- list(s = s_i, w = w_i)
  mcse <- list(s = s_i, w = w_i)
  for (m in metrics) for (o in outcomes) for (mth in methods) {
    st <- msd(df[[paste(m, o, mth, sep = "_")]])
    perf[[paste(m, o, mth, sep = "_")]] <- fmt(st["mean"], st["sd"])
    mcse[[paste(m, o, mth, sep = "_")]] <- round(unname(st["mcse"]), 4)
  }
  perf_rows[[sc]] <- perf
  mcse_rows[[sc]] <- mcse

  ## 3b) Density + coverage (JM all 200; LM bootstrap subset)
  cl_t <- msd(df$cov_trans_lm)
  cov_rows[[sc]] <- list(
    s = s_i, w = w_i,
    subjects_at_risk = round(mean(num(df$n_at_risk),     na.rm = TRUE), 1),
    events_trans     = round(mean(num(df$n_events_trans),na.rm = TRUE), 1),
    events_death     = round(mean(num(df$n_events_death),na.rm = TRUE), 1),
    cov_trans_JM = round(unname(msd(df$cov_trans_jm)["mean"]), 3),
    cov_death_JM = round(unname(msd(df$cov_death_jm)["mean"]), 3),
    cov_trans_LM = round(unname(cl_t["mean"]), 3),
    cov_death_LM = round(unname(msd(df$cov_death_lm)["mean"]), 3),
    n_reps_LM    = unname(cl_t["n"])     # number of replications LM coverage is based on
  )

  ## 3c) Mechanism diagnosis
  bt_jm <- num(df$bias_trans_jm); bd_jm <- num(df$bias_death_jm)
  bt_lm <- num(df$bias_trans_lm); bd_lm <- num(df$bias_death_lm)
  brk        <- which(bt_jm > break_thr)                 # which() drops NAs
  valid_idx  <- which(!is.na(bt_jm))
  stable_idx <- setdiff(valid_idx, brk)                  # safe against the empty-brk trap
  ev_t       <- mean(num(df$n_events_trans), na.rm = TRUE)
  nrep       <- length(valid_idx)
  mech_rows[[sc]] <- list(
    s = s_i, w = w_i,
    events_trans = round(ev_t, 1),
    EPV          = round(ev_t / n_params_surv, 2),
    r_JM = round(cor(bt_jm, bd_jm, use = "complete.obs"), 3),
    r_LM = round(cor(bt_lm, bd_lm, use = "complete.obs"), 3),
    pct_breaking      = round(100 * length(brk) / nrep, 1),
    break_bias_trans  = round(mean(bt_jm[brk]), 3),
    break_bias_death  = round(mean(bd_jm[brk]), 3),
    stable_bias_trans = round(mean(bt_jm[stable_idx]), 3),
    stable_bias_death = round(mean(bd_jm[stable_idx]), 3),
    converged_pct = round(100 * sum(num(df$converged), na.rm = TRUE) / nrep, 1),
    mean_Rhat     = round(mean(num(df$max_rhat), na.rm = TRUE), 3)
  )
}

to_df <- function(rows)
  do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
perf_df <- to_df(perf_rows)
mcse_df <- to_df(mcse_rows)
cov_df  <- to_df(cov_rows)
mech_df <- to_df(mech_rows)

## ---- 4) HUMAN-READABLE CONSOLE SUMMARY --------------------------------------
options(width = 220)
line <- function(t) cat("\n", strrep("=", 22), " ", t, " ", strrep("=", 22), "\n", sep = "")
line("PERFORMANCE  —  mean (SD), 200 replications  (Table 2)")
print(perf_df, row.names = FALSE)
line("COVERAGE & DATA DENSITY  (Table S3)")
print(cov_df, row.names = FALSE)
line("MECHANISM DIAGNOSIS  (Table S4)")
print(mech_df, row.names = FALSE)
line("MONTE CARLO SE — MCSE  (Table S5)")
print(mcse_df, row.names = FALSE)

## ---- 5) WRITE PUBLICATION-READY TABLES ---------------------------------
write.csv(perf_df, "OUT_table2_performance.csv",        row.names = FALSE)
write.csv(cov_df,  "OUT_tableS3_coverage_density.csv",  row.names = FALSE)
write.csv(mech_df, "OUT_tableS4_mechanism.csv",         row.names = FALSE)
write.csv(mcse_df, "OUT_tableS5_mcse.csv",              row.names = FALSE)
cat("\n[OK] Written: OUT_table2_performance.csv, OUT_tableS3_coverage_density.csv,",
    "OUT_tableS4_mechanism.csv, OUT_tableS5_mcse.csv\n")

## ---- EXPECTED VALUES (verified with Python; compare against your own output) ----
#  r_JM:  s2.5w2 = -0.937 | s2.5w4 = -0.962 | s5.5w2 = -0.975 | s5.5w4 = -0.984   (r_LM ~ 0)
#  %breaking: 7.5 / 15.5 / 19.0 / 19.0    | breaking bias (s5.5w4): trans +0.231, death -0.227
#  events_trans: 12.9 / 28.4 / 16.4 / 33.6    | converged: 58.5 / 58.5 / 62.0 / 62.0
#  cov_trans_JM ~0.75-0.79 ; cov_trans_LM ~0.92-0.95 ; cov_death_LM at s*_w4 ~0.63
