# Joint vs. Landmark – Competing Risks Simulation

## Contents
- **run_sim.R** – Main simulation script (compares Joint Modeling and Landmarking).
- **survfitLMCR.R** – Landmark prediction function (from Ferrer et al., 2019).
- **results/** – Output folder where simulation results will be saved.
- **LICENSE** – MIT license.

## How to Run
1. Extract all files into a single folder (e.g., `joint-vs-landmark-comprisk-sim/`).
2. Open R or RStudio and set the working directory to this folder:
   ```r
   setwd("C:/Users/username/joint-vs-landmark-comprisk-sim")
   source("run_sim.R")
3. After completion, the summary table will be saved as:
results/mc_summary_metrics.csv

## Requirements
### R version: 4.4.1 (tested under Windows 11 x64)
### Key R packages (tested versions):
 - JMbayes2 0.5-2
 - nlme 3.1-168
 - lme4 1.1-37
 - survival 3.8-3
 - riskRegression 2023.12.21
 - prodlim 2024.06.25
 - MASS 7.3-60.2
 - dplyr 1.1.4
 - timeROC 0.4
 - haven 2.5.4
Additional packages are automatically loaded as dependencies (Matrix, GLMMadaptive, etc.).

## References

van Houwelingen, H. C. (2007). Dynamic prediction by landmarking in event history analysis. Scandinavian Journal of Statistics, 34(1), 70–85. https://doi.org/10.1111/j.1467-9469.2006.00529.x
Rizopoulos, D., Miranda Afonso, P., & Papageorgiou, G. (2025). JMbayes2: Extended Joint Models for Longitudinal and Time-to-Event Data. R Package  Version 0.5-2. https://cran.r-project.org/package=JMbayes2
Ferrer, L., Putter, H., & Proust-Lima, C. (2019). Individual dynamic predictions using landmarking and joint modelling: Validation of estimators and robustness assessment. Statistical Methods in Medical Research, 28(12), 3649–3666. https://doi.org/10.1177/0962280218811837
