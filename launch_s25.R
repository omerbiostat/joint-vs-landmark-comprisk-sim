# =============================================================================
# launch_s25.R
# -----------------------------------------------------------------------------
# Runs the s = 2.5 scenarios (horizons w = 2 and w = 4 together).
# Start it in the FIRST RStudio session with source('launch_s25.R').
#
# Usage:
#   1. Open this file in RStudio
#   2. Paste into the console:
#        setwd("path/to/joint-vs-landmark-comprisk-sim")
#        source("launch_s25.R")
#   3. Wait several days. Keep the computer on.
#   4. If it crashes, re-run the same command; it resumes from the checkpoint.
# =============================================================================

cat("\n=== LAUNCH s=2.5 ===\n")
cat("Start: ", as.character(Sys.time()), "\n\n")

# Load the core functions
source("run_sim_v6_core.R")

# Run
results_s25 <- run_scenario_v6(
  s            = 2.5,
  w_list       = c(2, 4),
  M_total      = 200,
  CI_SUBSET    = 50,
  Mboot        = 200,
  n            = 1000,
  CensorTime   = 15,
  seed         = 2026,
  output_dir   = "results/s25",
  checkpoint_every = 20
)

cat("\nEnd: ", as.character(Sys.time()), "\n")
cat("\nOutput files:\n")
cat("  results/s25/results_s2.5_w2.csv\n")
cat("  results/s25/results_s2.5_w4.csv\n")
cat("  results/s25/checkpoint.rds  (complete)\n")
