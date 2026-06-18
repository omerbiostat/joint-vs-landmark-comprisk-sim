# Joint Modeling vs. Landmarking: A Simulation Study of Performance in Dynamic Prediction under Competing Risks

This repository contains the R code for a simulation study comparing **Joint Modeling (JM)** and **Landmarking (LM)** for dynamic cumulative-incidence prediction under competing risks (transplantation and death), with the data-generating parameters calibrated to the PBC2 cohort. It accompanies the manuscript of the same title (under review at *PLOS ONE*).

## Contents

- **run_sim_v6_core.R** — Core simulation functions. Generates one replication of longitudinal and competing-risks survival data, fits both the joint model and the landmark model, computes the performance metrics, and runs the full 200-replication scenario with checkpointing. *Not run directly* — it is loaded by the launcher scripts.
- **launch_s25.R** — Launcher for the early-landmark scenarios (s = 2.5; horizons w = 2 and 4). Writes its output to `results/s25/`.
- **launch_s55.R** — Launcher for the late-landmark scenarios (s = 5.5; horizons w = 2 and 4). Writes its output to `results/s55/`.
- **survfitLMCR.R** — Landmark cause-specific prediction function with bootstrap confidence intervals, based on the two-stage approach of Ferrer et al. (2019).
- **summary_statistics.R** — Post-processing script (base R) that reads the four result CSVs and writes the publication-ready summary tables, including Monte Carlo standard errors.
- **figures.R** — Post-processing script that produces the manuscript figures from the four result CSVs.
- **results/** — Output folder where the per-scenario result files are written.
- **LICENSE** — MIT license.

## How to Run

1. Extract all files into a single folder (e.g., `joint-vs-landmark-comprisk-sim/`).
2. The two landmark times can be run in parallel in **two separate R/RStudio sessions**. In each session, set the working directory to this folder:
   ```r
   setwd("path/to/joint-vs-landmark-comprisk-sim")
   ```
3. In the first session, run the early-landmark scenarios:
   ```r
   source("launch_s25.R")
   ```
   In the second session, run the late-landmark scenarios:
   ```r
   source("launch_s55.R")
   ```
   (Each launcher loads `run_sim_v6_core.R`, which in turn sources `survfitLMCR.R`.)
4. **This is computationally intensive.** Each landmark runs 200 Monte Carlo replications at n = 1000; on a typical desktop this takes on the order of several days per landmark, driven by the Bayesian joint-model fitting.
5. **Checkpointing / restart.** Progress is saved every 20 replications (`results/s25/checkpoint.rds`, `results/s55/checkpoint.rds`). If a session is interrupted, simply re-run the same `source(...)` command; it resumes from the last checkpoint.
6. On completion, the per-scenario results are written to:
   ```
   results/s25/results_s2.5_w2.csv
   results/s25/results_s2.5_w4.csv
   results/s55/results_s5.5_w2.csv
   results/s55/results_s5.5_w4.csv
   ```
   Each row is one replication. Columns include bias, RMSE, Brier score, and time-dependent AUC for transplantation and death under both JM and LM, 95% interval coverage, the average risk-set size and event counts, and convergence diagnostics (maximum R-hat and a convergence flag).

**Reproducibility.** Deterministic seeds are used (seed = 2026 for s = 2.5 and seed = 5555 for s = 5.5), with a distinct seed per replication, so the results are exactly reproducible.

## Post-processing

Once the four result CSVs have been produced, the two post-processing scripts reproduce the tables and figures in the manuscript. Before running either one, open it and set `data_dir` (near the top) to the folder that holds the four CSV files.

- **summary_statistics.R** (base R only) — Summarises the 200 replications and writes four publication-ready tables: `OUT_table2_performance.csv` (Table 2; mean and SD per metric for JM and LM), `OUT_tableS3_coverage_density.csv` (S3 Table; data density and interval coverage), `OUT_tableS4_mechanism.csv` (S4 Table; mechanism-diagnosis and convergence metrics), and `OUT_tableS5_mcse.csv` (S5 Table; Monte Carlo standard errors). It also prints a readable summary to the console.
- **figures.R** — Produces the manuscript figures as 300-dpi TIFF and PNG: `Figure1_performance`, `Figure2_mechanism`, `FigureS1_coverage`, and `FigureS2_bias_distributions`. This script additionally requires `ggplot2`, `dplyr`, `tidyr`, and `patchwork`.

Run them with:
```r
source("summary_statistics.R")
source("figures.R")
```

## Requirements

**R version:** 4.4.1

**Key R packages:**
- JMbayes2 0.5-2
- nlme 3.1-168
- survival 3.8-3
- riskRegression 2023.12.21
- prodlim 2024.06.25
- MASS 7.3-60.2
- dplyr 1.1.4
- splines (part of base R)

## Notes on this version

This version corresponds to the revised manuscript and reflects the following choices: a single, clinically realistic cohort size (n = 1000); prediction horizons of 2 and 4 years; a landmark predictor based on the current biomarker level (no slope term); and Landmarking uncertainty obtained by bootstrap. It supersedes the earlier single `run_sim.R` script, which has been replaced by the core/launcher structure above.

## References

- Ferrer L, Putter H, Proust-Lima C. Individual dynamic predictions using landmarking and joint modelling: Validation of estimators and robustness assessment. *Statistical Methods in Medical Research.* 2019;28(12):3649–3666. *(reference for the landmark prediction function)*
