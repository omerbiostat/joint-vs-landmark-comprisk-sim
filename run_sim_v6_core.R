# =============================================================================
# run_sim_v6_core.R
# -----------------------------------------------------------------------------
# Core functions for the v6 simulation.
# Do NOT run this file directly - the launcher files (launch_s25.R,
# launch_s55.R) load it via source().
#
# Provides:
#   simulate_joint_data()       : generates data for one replication
#   simulate_once_v6()          : computes one replication (JM + LM + scores)
#   run_scenario_v6()           : runs 200 replications with checkpointing
#
# Differences from v5:
#   - n = 1000 (fixed)
#   - JM settings: 25000 / 5000 / 5 / 3 chains
#   - no slope in LM (Reviewer #1)
#   - JM-shared structure: a single JM fit, separate prediction for each horizon
#   - LM-shared structure: a single LME fit, separate Cox + bootstrap for each horizon
#   - IPCW landmark-shifted Score (following Blanche et al. 2015)
#   - Checkpoint system: save an RDS every 20 replications
#   - Restart-aware: an interrupted run resumes from where it stopped
# =============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(JMbayes2)
  library(nlme)
  library(MASS)
  library(dplyr)
  library(splines)
  library(riskRegression)
  library(prodlim)
})

source("survfitLMCR.R")

# =============================================================================
# DGP parameters and baseline hazard objects - GLOBAL
# =============================================================================
.DGP <- new.env()

with(.DGP, {
  truebeta  <- c(`(Intercept)` = 10.63765,
                 year          = 0.2629633,
                 drug          = -0.09780603,
                 `year:drug`   = -0.01708882)
  sigma_e   <- 1.058143
  Dmat      <- matrix(c(0.6152592, 0.007199312, 0.007199312, 0.09915757), 2, 2)
  truegamma <- c(transplanted = -0.4647676, dead = 0.3482824)
  truealpha <- c(transplanted = 0.1553865,  dead = 0.2506399)

  T_max         <- 15
  n_knots       <- 5
  knots         <- seq(0, T_max, length = n_knots + 2)[-c(1, n_knots + 2)]
  K             <- length(knots) + 3 + 1
  gamma0        <- seq(-3, 0, length.out = K)
  n_intervals   <- 1000
  cuts          <- seq(0, T_max, length.out = n_intervals + 1)
  delta_t       <- diff(cuts)
  midpt         <- (cuts[-1] + cuts[-length(cuts)]) / 2
  f_tr          <- 0.014
  f_de          <- 0.0132

  bs_midpt  <- bs(midpt, knots = knots, degree = 3, intercept = TRUE)
  h_base_tr <- as.numeric(f_tr * exp(bs_midpt %*% gamma0))
  h_base_de <- as.numeric(f_de * exp(bs_midpt %*% gamma0))
})


# =============================================================================
# Generate a single event time
# =============================================================================
gen_event_time <- function(a_i, d_i, drug_i, gamma_k, alpha_k, h_base_cause) {
  a_i      <- as.numeric(a_i)
  d_i      <- as.numeric(d_i)
  drug_i   <- as.numeric(drug_i)
  gamma_k  <- as.numeric(gamma_k)
  alpha_k  <- as.numeric(alpha_k)

  U <- as.numeric(runif(1))
  H_target <- -log(U)

  lin_pred <- gamma_k * drug_i + alpha_k * (a_i + d_i * .DGP$midpt)
  lin_pred <- pmin(lin_pred, 20)
  eff_haz  <- h_base_cause * exp(lin_pred)
  H_cum    <- cumsum(eff_haz * .DGP$delta_t)

  j <- which(H_cum >= H_target)[1]
  if (is.na(j)) return(.DGP$T_max + 10)

  H_prev <- if (j > 1) H_cum[j - 1] else 0
  t_prev <- .DGP$cuts[j]
  t_prev + (H_target - H_prev) / eff_haz[j]
}


# =============================================================================
# Data generation (for n patients)
# =============================================================================
simulate_joint_data <- function(n, CensorTime) {

  drug_vec <- rbinom(n, 1, 0.5)
  b_mat    <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = .DGP$Dmat)

  long_list  <- vector("list", n)
  surv_list  <- vector("list", n)
  truth_list <- vector("list", n)

  for (i in seq_len(n)) {
    a_i <- .DGP$truebeta[["(Intercept)"]] + b_mat[i, 1] +
           drug_vec[i] * .DGP$truebeta[["drug"]]
    d_i <- .DGP$truebeta[["year"]] + b_mat[i, 2] +
           drug_vec[i] * .DGP$truebeta[["year:drug"]]
    drug_i <- drug_vec[i]

    T1 <- gen_event_time(a_i, d_i, drug_i,
                         gamma_k = .DGP$truegamma[["transplanted"]],
                         alpha_k = .DGP$truealpha[["transplanted"]],
                         h_base_cause = .DGP$h_base_tr)
    T2 <- gen_event_time(a_i, d_i, drug_i,
                         gamma_k = .DGP$truegamma[["dead"]],
                         alpha_k = .DGP$truealpha[["dead"]],
                         h_base_cause = .DGP$h_base_de)

    T_true_event <- min(T1, T2)
    true_status  <- if (T1 < T2) 1 else 2

    Cens_admin  <- CensorTime
    Cens_random <- rexp(1, rate = 0.05)
    Cens_i      <- min(Cens_admin, Cens_random)
    T_i     <- min(T_true_event, Cens_i)
    event_i <- if (Cens_i < T_true_event) 0 else true_status

    protocol_visits <- seq(0, 15, by = 0.8)
    possible_visits <- protocol_visits[protocol_visits <= T_i]
    if (length(possible_visits) == 0) possible_visits <- 0

    if (length(possible_visits) > 1) {
      kept_indices  <- rbinom(length(possible_visits) - 1, 1, 0.85)
      actual_visits <- c(0, possible_visits[-1][kept_indices == 1])
    } else {
      actual_visits <- possible_visits
    }

    t_ij <- actual_visits + rnorm(length(actual_visits), 0, 0.1)
    if (length(t_ij) > 0) t_ij[1] <- 0
    t_ij <- sort(t_ij)
    t_ij <- t_ij[t_ij <= T_i]
    t_ij <- pmax(0, t_ij)

    mu_ij <- a_i + d_i * t_ij
    Y_ij  <- mu_ij + rnorm(length(t_ij), 0, .DGP$sigma_e)

    long_list[[i]]  <- data.frame(id = i, years = T_i, time = t_ij,
                                  Y = Y_ij, drug = drug_i)
    surv_list[[i]]  <- data.frame(id = i, years = T_i, time = min(t_ij),
                                  Y = Y_ij[1], event = event_i, drug = drug_i)
    truth_list[[i]] <- data.frame(id = i, a_i = a_i, d_i = d_i)
  }

  list(long  = do.call(rbind, long_list),
       surv  = do.call(rbind, surv_list),
       truth = do.call(rbind, truth_list))
}


# =============================================================================
# True conditional CIF (analytic)
# =============================================================================
true_cif_cond <- function(a_i, d_i, drug_i, s, t_star,
                          gamma_tr, alpha_tr, gamma_de, alpha_de) {

  midpt   <- .DGP$midpt
  delta_t <- .DGP$delta_t
  cuts    <- .DGP$cuts
  h_tr    <- .DGP$h_base_tr
  h_de    <- .DGP$h_base_de

  lin1 <- pmin(gamma_tr * drug_i + alpha_tr * (a_i + d_i * midpt), 20)
  lin2 <- pmin(gamma_de * drug_i + alpha_de * (a_i + d_i * midpt), 20)
  h1   <- h_tr * exp(lin1)
  h2   <- h_de * exp(lin2)
  hs   <- h1 + h2

  eval_at <- function(t) {
    if (t <= 0)        return(list(S = 1, CIF1 = 0, CIF2 = 0))
    if (t >= max(cuts)) t <- max(cuts)
    j  <- findInterval(t, cuts, rightmost.closed = TRUE)
    j  <- min(j, length(delta_t))
    dt <- t - cuts[j]

    if (j == 1) {
      S_before    <- 1
      CIF1_before <- 0
      CIF2_before <- 0
    } else {
      H_inc   <- hs * delta_t
      H_end   <- cumsum(H_inc)
      S_start <- c(1, exp(-H_end[-length(H_end)]))
      inc1    <- S_start * (h1/hs) * (1 - exp(-hs * delta_t))
      inc2    <- S_start * (h2/hs) * (1 - exp(-hs * delta_t))
      S_before    <- exp(-H_end[j - 1])
      CIF1_before <- cumsum(inc1)[j - 1]
      CIF2_before <- cumsum(inc2)[j - 1]
    }

    hz <- hs[j]
    if (hz <= 0) {
      list(S = S_before, CIF1 = CIF1_before, CIF2 = CIF2_before)
    } else {
      fac    <- 1 - exp(-hz * dt)
      CIF1_t <- CIF1_before + S_before * (h1[j]/hz) * fac
      CIF2_t <- CIF2_before + S_before * (h2[j]/hz) * fac
      list(S = exp(-(-log(S_before) + hz * dt)),
           CIF1 = CIF1_t, CIF2 = CIF2_t)
    }
  }

  at_s  <- eval_at(s)
  at_ts <- eval_at(t_star)

  S_s <- at_s$S
  if (S_s <= 0) return(c(CIF1 = NA_real_, CIF2 = NA_real_))

  CIF1_c <- max(0, min(1, (at_ts$CIF1 - at_s$CIF1) / S_s))
  CIF2_c <- max(0, min(1, (at_ts$CIF2 - at_s$CIF2) / S_s))
  tot <- CIF1_c + CIF2_c
  if (tot > 1) { CIF1_c <- CIF1_c/tot; CIF2_c <- CIF2_c/tot }
  c(CIF1 = CIF1_c, CIF2 = CIF2_c)
}


# =============================================================================
# LM model: with parameters tLM and thor
# No slope from the Cox model (only X + level)
# =============================================================================
do_LM <- function(sim.long, sim.surv, tLM, thor,
                  simulate_flag = FALSE, M = 200) {
  tpred  <- tLM + thor
  Ri_LM  <- sim.surv$id[sim.surv$years > tLM]
  nLM    <- length(Ri_LM)
  LMlong <- subset(sim.long, time < tLM & id %in% Ri_LM)
  LMlong <- LMlong %>% rename(X = drug)
  LMsurv <- subset(sim.surv, id %in% Ri_LM)
  LMsurv <- LMsurv %>% rename(X = drug)
  LMsurv$Rec   <- as.integer(LMsurv$event == 1)
  LMsurv$Death <- as.integer(LMsurv$event == 2)

  lmeFit <- try(lme(fixed = Y ~ time * X, data = LMlong,
                    random = ~ time | id,
                    control = lmeControl(opt = "optim",
                                         msMaxIter = 1000, msMaxEval = 1000)),
                silent = TRUE)
  if (inherits(lmeFit, "try-error")) {
    lmeFit <- lme(fixed = Y ~ time * X, data = LMlong, random = ~ 1 | id,
                  control = lmeControl(opt = "optim"))
  }

  b       <- ranef(lmeFit)
  sigma_l <- lmeFit$sigma
  Dmat_l  <- lapply(pdMatrix(lmeFit$modelStruct$reStruct), "*", sigma_l^2)[[1]]
  betas   <- fixef(lmeFit)

  Data.s <- data.frame(id = unique(LMlong$id), time = tLM,
                       X = LMlong[!duplicated(LMlong$id), "X"])
  Xlong.s <- model.matrix(~ time * X, Data.s)
  Z.s     <- model.matrix(~ time, Data.s)
  LMsurv$level <- as.vector(c(Xlong.s %*% betas) + rowSums(Z.s * b))

  LMsurv$tRecAC   <- pmin(LMsurv$years, tpred)
  LMsurv$RecAC    <- LMsurv$Rec
  LMsurv$RecAC[LMsurv$years > tpred] <- 0
  LMsurv$tDeathAC <- pmin(LMsurv$years, tpred)
  LMsurv$DeathAC  <- LMsurv$Death
  LMsurv$DeathAC[LMsurv$years > tpred] <- 0

  coxFit1 <- coxph(Surv(tRecAC,   RecAC)   ~ X + level,
                   data = LMsurv, x = TRUE, y = TRUE)
  coxFit2 <- coxph(Surv(tDeathAC, DeathAC) ~ X + level,
                   data = LMsurv, x = TRUE, y = TRUE)

  Xsurv_pred <- model.matrix(~ 0 + X + level, LMsurv)
  HR1 <- as.numeric(exp(Xsurv_pred %*% coxFit1$coef))
  HR2 <- as.numeric(exp(Xsurv_pred %*% coxFit2$coef))
  bh1 <- basehaz(coxFit1, centered = FALSE)
  bh2 <- basehaz(coxFit2, centered = FALSE)

  toci <- function(bh1, bh2, HR1, HR2, tpred){
    h1 <- bh1; names(h1)[1] <- "Haz"; h1$Haz <- h1$Haz * HR1; h1$cause <- 1
    h2 <- bh2; names(h2)[1] <- "Haz"; h2$Haz <- h2$Haz * HR2; h2$cause <- 2
    Haz <- rbind(h1, h2)
    CI  <- CumInc(Haz)
    idx <- sum(CI$time <= tpred)
    return(CI[idx, ])
  }

  ci <- matrix(NA, nLM, 2)
  for (i in seq_len(nLM)) {
    ci[i, ] <- as.numeric(toci(bh1, bh2, HR1[i], HR2[i], tpred))[2:3]
  }

  res.MC <- list(Rec = NULL, Death = NULL)

  if (simulate_flag) {
    tt     <- matrix(coxFit1$y, ncol = 2)[, 1]
    Delta1 <- matrix(coxFit1$y, ncol = 2)[, 2]
    Delta2 <- matrix(coxFit2$y, ncol = 2)[, 2]
    aV <- lmeFit$apVar
    if (!is.null(aV) && !(is.character(aV) && length(aV) == 1) &&
        all(is.finite(aV))) {
      Pars   <- attr(aV, "Pars")
      varFix <- lmeFit$varFix
      nbetas <- length(betas)
      nP     <- length(Pars)
      mat <- matrix(0, nbetas + nP, nbetas + nP)
      mat[seq_len(nbetas), seq_len(nbetas)]           <- varFix
      mat[nbetas + seq_len(nP), nbetas + seq_len(nP)] <- aV

      coef_cox1.MC <- mvrnorm(M, coxFit1$coef, coxFit1$var)
      coef_cox2.MC <- mvrnorm(M, coxFit2$coef, coxFit2$var)
      coef_long.MC <- mvrnorm(M, c(betas, Pars), mat)

      ci.MC <- replicate(M, matrix(NA, nrow = nLM, ncol = 2), simplify = FALSE)
      Xlong_l <- split(data.frame(model.matrix(~ time * X, LMlong)), LMlong$id)
      Z_l     <- split(data.frame(model.matrix(~ time, LMlong)),     LMlong$id)
      Y_l     <- split(LMlong$Y, LMlong$id)

      for (l in seq_len(M)) {
        betas.MC <- coef_long.MC[l, seq_len(nbetas)]
        sigma.MC <- exp(coef_long.MC[l, nbetas + nP])
        lmeSt <- lmeFit$modelStruct
        lmeSt$reStruct[[1]] <- pdNatural(lmeSt$reStruct[[1]])
        coef(lmeSt) <- coef_long.MC[l, -c(seq_len(nbetas), nbetas + nP)]
        Pars_D.MC   <- coef(lmeSt, unconstrained = FALSE)
        D.MC <- matrix(c(Pars_D.MC[1], rep(Pars_D.MC[3], 2), Pars_D.MC[2]), 2, 2)
        diag(D.MC) <- diag(D.MC)^2
        D.MC[upper.tri(D.MC, diag = FALSE)] <-
          D.MC[lower.tri(D.MC, diag = FALSE)] <- Reduce("*", Pars_D.MC)

        V.MC <- lapply(Z_l, function(x){
                  as.matrix(x) %*% D.MC %*% t(as.matrix(x)) +
                  diag(sigma.MC^2, nrow(as.matrix(x)))})
        b.MC <- lapply(seq_along(Xlong_l), function(i){
                  D.MC %*% t(as.matrix(Z_l[[i]])) %*% solve(V.MC[[i]]) %*%
                  (Y_l[[i]] - as.matrix(Xlong_l[[i]]) %*% betas.MC)})
        b.MC <- matrix(unlist(b.MC), ncol = 2, byrow = TRUE)

        LMsurv.MC <- LMsurv
        LMsurv.MC$level <- as.vector(c(Xlong.s %*% betas.MC) +
                                     rowSums(Z.s * b.MC))
        Xsurv.MC <- model.matrix(~ 0 + X + level, LMsurv.MC)
        HR1.MC <- as.numeric(exp(Xsurv.MC %*% coef_cox1.MC[l, ]))
        HR2.MC <- as.numeric(exp(Xsurv.MC %*% coef_cox2.MC[l, ]))

        set.seed(l)
        nu_event <- 4 * rbeta(nLM, 1/2, 3/2)
        haz1.PR <- haz2.PR <- N1.PR <- N2.PR <- numeric(nLM)
        for (i in 1:nLM) {
          haz1.PR[i] <- mean(nu_event * as.numeric(tt >= tt[i]) *
                             exp(c(Xsurv.MC %*% coef_cox1.MC[l, ])))
          haz2.PR[i] <- mean(nu_event * as.numeric(tt >= tt[i]) *
                             exp(c(Xsurv.MC %*% coef_cox2.MC[l, ])))
          N1.PR[i]   <- mean(nu_event * as.numeric(tt <= tt[i]) * Delta1)
          N2.PR[i]   <- mean(nu_event * as.numeric(tt <= tt[i]) * Delta2)
        }
        ord <- order(tt)
        Haz1.PR <- data.frame(hazard = cumsum(diff(c(0, N1.PR[ord])) /
                                              haz1.PR[ord]),
                              time   = tt[ord])
        Haz2.PR <- data.frame(hazard = cumsum(diff(c(0, N2.PR[ord])) /
                                              haz2.PR[ord]),
                              time   = tt[ord])

        for (i in 1:nLM) {
          ci.MC[[l]][i, ] <- as.numeric(toci(Haz1.PR, Haz2.PR,
                                             HR1.MC[i], HR2.MC[i], tpred))[2:3]
        }
      }

      res.MC <- replicate(2, matrix(NA, nrow = nLM, ncol = 6), simplify = FALSE)
      for (k in 1:2) {
        cik <- sapply(ci.MC, function(x) x[, k])
        if (is.null(dim(cik))) cik <- matrix(cik, nrow = nLM)
        for (i in seq_len(nLM)) {
          vals <- cik[i, ]
          res.MC[[k]][i, ] <- c(tLM, thor, mean(vals, na.rm = TRUE),
                                median(vals, na.rm = TRUE),
                                quantile(vals, 0.025, na.rm = TRUE),
                                quantile(vals, 0.975, na.rm = TRUE))
        }
        colnames(res.MC[[k]]) <- c("tLM", "thor", "Mean", "Median",
                                   "Lower (2.5%)", "Upper (97.5%)")
        rownames(res.MC[[k]]) <- paste0("ID", LMsurv$id)
      }
      names(res.MC) <- c("Rec", "Death")
    }
  }

  list(id_vec = LMsurv$id,
       trans  = ci[, 1],
       death  = ci[, 2],
       res.MC = res.MC)
}


# =============================================================================
# Single replication: one landmark s, two horizons (w_list)
# JM is fitted once, separate prediction for each horizon
# LM likewise (LME once, Cox twice)
# =============================================================================
simulate_once_v6 <- function(s, w_list,
                             n = 1000, CensorTime = 15,
                             compute_lm_ci = FALSE,
                             Mboot = 200) {

  # --- Generate data ---
  sim <- simulate_joint_data(n, CensorTime)
  sim.long  <- sim$long
  sim.surv  <- sim$surv
  sim.truth <- sim$truth

  # --- JM fit (ONCE, for landmark s) ---
  simCR <- crisk_setup(data = sim.surv, statusVar = "event",
                       censLevel = 0, nameStrata = "CR")

  CoxFit_CR <- coxph(Surv(years, status2) ~ drug * strata(CR),
                     data = simCR, x = TRUE, y = TRUE)
  fm_Y <- lme(Y ~ time * drug, random = ~ time | id, data = sim.long,
              control = lmeControl(opt = "optim"))

  jFit_CR <- jm(Surv_object   = CoxFit_CR,
                Mixed_objects = list(fm_Y),
                time_var      = "time",
                functional_forms = list("Y" = ~ value(Y) : CR),
                id_var        = "id",
                n_iter   = 25000L,
                n_burnin = 5000L,
                n_thin   = 5L,
                n_chains = 3L,
                cores    = 3)

  rhat_vals <- unlist(jFit_CR$statistics$Rhat)
  rhat_vals <- rhat_vals[is.finite(rhat_vals)]
  max_rhat  <- if (length(rhat_vals)) max(rhat_vals) else NA_real_
  converged <- isTRUE(max_rhat < 1.1)

  Ri <- unique(sim.surv$id[sim.surv$years > s])

  # --- Metrics for each horizon ---
  results_list <- list()

  for (w in w_list) {
    t_star <- s + w

    # JM predictions
    jm_df_w <- data.frame(id = Ri, joint_trans = NA, trans_low = NA, trans_upp = NA,
                          joint_death = NA, death_low = NA, death_upp = NA)
    for (j in seq_along(Ri)) {
      idj <- Ri[j]
      ND_long  <- sim.long[sim.long$id == idj & sim.long$time < s, ]
      ND_event <- simCR[simCR$id == idj, ]
      ND_event$status2 <- 0
      ND_event$years   <- s
      ND <- list(newdataL = ND_long, newdataE = ND_event)

      P <- predict(jFit_CR, newdata = ND, return_newdata = TRUE,
                   times = t_star, process = "event")

      jm_df_w$joint_trans[j] <- with(P, P$pred_CIF[`_strata` == 2 & years == t_star])
      jm_df_w$trans_low[j]   <- with(P, P$low_CIF[`_strata`  == 2 & years == t_star])
      jm_df_w$trans_upp[j]   <- with(P, P$upp_CIF[`_strata`  == 2 & years == t_star])
      jm_df_w$joint_death[j] <- with(P, P$pred_CIF[`_strata` == 1 & years == t_star])
      jm_df_w$death_low[j]   <- with(P, P$low_CIF[`_strata`  == 1 & years == t_star])
      jm_df_w$death_upp[j]   <- with(P, P$upp_CIF[`_strata`  == 1 & years == t_star])
    }

    # LM point predictions
    lm_det <- do_LM(sim.long, sim.surv, s, w, simulate_flag = FALSE)
    lm_df  <- data.frame(id = lm_det$id_vec,
                         Rec_Mean = lm_det$trans,
                         Death_Mean = lm_det$death)

    out <- merge(jm_df_w, lm_df, by = "id")

    # True CIF
    truth_sub <- merge(sim.truth, sim.surv[, c("id", "drug")], by = "id", sort = FALSE)
    truth_sub <- truth_sub[match(out$id, truth_sub$id), ]

    true_mat <- mapply(
      FUN = function(a, d, dr) true_cif_cond(
        a_i = a, d_i = d, drug_i = dr, s = s, t_star = t_star,
        gamma_tr = .DGP$truegamma[["transplanted"]],
        alpha_tr = .DGP$truealpha[["transplanted"]],
        gamma_de = .DGP$truegamma[["dead"]],
        alpha_de = .DGP$truealpha[["dead"]]),
      a = truth_sub$a_i, d = truth_sub$d_i, dr = truth_sub$drug)

    true_trans <- as.numeric(true_mat["CIF1", ])
    true_death <- as.numeric(true_mat["CIF2", ])

    # Bias and RMSE
    bias_trans_jm <- mean(out$joint_trans - true_trans, na.rm = TRUE)
    rmse_trans_jm <- sqrt(mean((out$joint_trans - true_trans)^2, na.rm = TRUE))
    bias_death_jm <- mean(out$joint_death - true_death, na.rm = TRUE)
    rmse_death_jm <- sqrt(mean((out$joint_death - true_death)^2, na.rm = TRUE))

    bias_trans_lm <- mean(out$Rec_Mean - true_trans, na.rm = TRUE)
    rmse_trans_lm <- sqrt(mean((out$Rec_Mean - true_trans)^2, na.rm = TRUE))
    bias_death_lm <- mean(out$Death_Mean - true_death, na.rm = TRUE)
    rmse_death_lm <- sqrt(mean((out$Death_Mean - true_death)^2, na.rm = TRUE))

    # JM coverage
    cov_trans_jm <- mean(true_trans >= out$trans_low &
                         true_trans <= out$trans_upp, na.rm = TRUE)
    cov_death_jm <- mean(true_death >= out$death_low &
                         true_death <= out$death_upp, na.rm = TRUE)

    # Score IPCW (landmark-shifted)
    evaldat <- sim.surv[match(out$id, sim.surv$id), c("years", "event")]
    evaldat$lm_time  <- pmin(evaldat$years - s, w)
    evaldat$lm_event <- evaldat$event
    evaldat$lm_event[evaldat$years - s > w] <- 0

    sc_tr_b <- Score(list(JM = out$joint_trans, LM = out$Rec_Mean),
                     formula = Hist(lm_time, lm_event) ~ 1, data = evaldat,
                     times = w, cause = 1, metrics = "brier",
                     cens.model = "km", conservative = TRUE)
    sc_de_b <- Score(list(JM = out$joint_death, LM = out$Death_Mean),
                     formula = Hist(lm_time, lm_event) ~ 1, data = evaldat,
                     times = w, cause = 2, metrics = "brier",
                     cens.model = "km", conservative = TRUE)
    sc_tr_a <- Score(list(JM = out$joint_trans, LM = out$Rec_Mean),
                     formula = Hist(lm_time, lm_event) ~ 1, data = evaldat,
                     times = w, cause = 1, metrics = "auc",
                     cens.model = "km", conservative = TRUE)
    sc_de_a <- Score(list(JM = out$joint_death, LM = out$Death_Mean),
                     formula = Hist(lm_time, lm_event) ~ 1, data = evaldat,
                     times = w, cause = 2, metrics = "auc",
                     cens.model = "km", conservative = TRUE)

    brier_trans_jm <- subset(sc_tr_b$Brier$score, model == "JM")$Brier[1]
    brier_trans_lm <- subset(sc_tr_b$Brier$score, model == "LM")$Brier[1]
    brier_death_jm <- subset(sc_de_b$Brier$score, model == "JM")$Brier[1]
    brier_death_lm <- subset(sc_de_b$Brier$score, model == "LM")$Brier[1]
    auc_trans_jm   <- subset(sc_tr_a$AUC$score,   model == "JM")$AUC[1]
    auc_trans_lm   <- subset(sc_tr_a$AUC$score,   model == "LM")$AUC[1]
    auc_death_jm   <- subset(sc_de_a$AUC$score,   model == "JM")$AUC[1]
    auc_death_lm   <- subset(sc_de_a$AUC$score,   model == "LM")$AUC[1]

    # LM coverage (only when compute_lm_ci=TRUE)
    cov_trans_lm <- NA_real_
    cov_death_lm <- NA_real_
    if (compute_lm_ci) {
      lm_ci_w <- try(do_LM(sim.long, sim.surv, s, w,
                           simulate_flag = TRUE, M = Mboot),
                     silent = TRUE)
      if (!inherits(lm_ci_w, "try-error") &&
          !is.null(lm_ci_w$res.MC$Rec) &&
          !is.null(lm_ci_w$res.MC$Death)) {
        ci_rec   <- lm_ci_w$res.MC$Rec
        ci_death <- lm_ci_w$res.MC$Death
        ci_ids   <- as.integer(gsub("ID *", "", rownames(ci_rec)))
        m_idx    <- match(out$id, ci_ids)
        cov_trans_lm <- mean(true_trans >= ci_rec[m_idx,   "Lower (2.5%)"] &
                             true_trans <= ci_rec[m_idx,   "Upper (97.5%)"], na.rm = TRUE)
        cov_death_lm <- mean(true_death >= ci_death[m_idx, "Lower (2.5%)"] &
                             true_death <= ci_death[m_idx, "Upper (97.5%)"], na.rm = TRUE)
      }
    }

    n_event_trans_window <- sum(sim.surv$years > s & sim.surv$years <= t_star &
                                sim.surv$event == 1)
    n_event_death_window <- sum(sim.surv$years > s & sim.surv$years <= t_star &
                                sim.surv$event == 2)

    results_list[[as.character(w)]] <- c(
      bias_trans_jm  = bias_trans_jm,
      rmse_trans_jm  = rmse_trans_jm,
      bias_death_jm  = bias_death_jm,
      rmse_death_jm  = rmse_death_jm,
      bias_trans_lm  = bias_trans_lm,
      rmse_trans_lm  = rmse_trans_lm,
      bias_death_lm  = bias_death_lm,
      rmse_death_lm  = rmse_death_lm,
      brier_trans_jm = brier_trans_jm,
      brier_trans_lm = brier_trans_lm,
      brier_death_jm = brier_death_jm,
      brier_death_lm = brier_death_lm,
      auc_trans_jm   = auc_trans_jm,
      auc_trans_lm   = auc_trans_lm,
      auc_death_jm   = auc_death_jm,
      auc_death_lm   = auc_death_lm,
      cov_trans_jm   = cov_trans_jm,
      cov_death_jm   = cov_death_jm,
      cov_trans_lm   = cov_trans_lm,
      cov_death_lm   = cov_death_lm,
      n_at_risk      = length(Ri),
      n_events_trans = n_event_trans_window,
      n_events_death = n_event_death_window,
      max_rhat       = max_rhat,
      converged      = as.numeric(converged)
    )
  }

  results_list
}


# =============================================================================
# Main scenario-running function (checkpoint system)
# =============================================================================
run_scenario_v6 <- function(s, w_list,
                            M_total = 200,
                            CI_SUBSET = 50,
                            Mboot = 200,
                            n = 1000,
                            CensorTime = 15,
                            seed = 2026,
                            output_dir,
                            checkpoint_every = 20) {

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  metric_names <- c(
    "bias_trans_jm", "rmse_trans_jm", "bias_death_jm", "rmse_death_jm",
    "bias_trans_lm", "rmse_trans_lm", "bias_death_lm", "rmse_death_lm",
    "brier_trans_jm", "brier_trans_lm", "brier_death_jm", "brier_death_lm",
    "auc_trans_jm", "auc_trans_lm", "auc_death_jm", "auc_death_lm",
    "cov_trans_jm", "cov_death_jm", "cov_trans_lm", "cov_death_lm",
    "n_at_risk", "n_events_trans", "n_events_death",
    "max_rhat", "converged"
  )

  # A separate matrix for each w
  results <- lapply(w_list, function(w) {
    matrix(NA_real_, nrow = M_total, ncol = length(metric_names),
           dimnames = list(NULL, metric_names))
  })
  names(results) <- as.character(w_list)

  # CI-subset replication indices (deterministic)
  set.seed(seed)
  ci_rep_idx <- sort(sample.int(M_total, size = min(CI_SUBSET, M_total)))

  # LOAD CHECKPOINT (if present)
  checkpoint_file <- file.path(output_dir, "checkpoint.rds")
  start_idx <- 1
  if (file.exists(checkpoint_file)) {
    cat("Checkpoint found, loading...\n")
    cp <- readRDS(checkpoint_file)
    results <- cp$results
    start_idx <- cp$next_idx
    cat(sprintf("Resuming from replication %d.\n", start_idx))
  }

  if (start_idx > M_total) {
    cat("This scenario is already complete. Exiting.\n")
    return(results)
  }

  cat(sprintf("\n==== Scenario: s=%g, w_list=(%s), M=%d ====\n",
              s, paste(w_list, collapse = ","), M_total))
  cat(sprintf("Output dir: %s\n", output_dir))
  cat(sprintf("CI subset replications: %d / %d\n", length(ci_rep_idx), M_total))
  cat(sprintf("Checkpoint every %d replications\n", checkpoint_every))
  cat(sprintf("Start index: %d\n\n", start_idx))

  total_start <- Sys.time()

  for (m in start_idx:M_total) {
    iter_start <- Sys.time()
    do_ci <- m %in% ci_rep_idx

    cat(sprintf("[%s] Replication %d / %d%s\n",
                format(Sys.time(), "%H:%M:%S"), m, M_total,
                if (do_ci) "  [+ LM CI]" else ""))

    # A separate seed per replication (reproducibility)
    set.seed(seed + m)

    res_m <- try(simulate_once_v6(s = s, w_list = w_list,
                                  n = n, CensorTime = CensorTime,
                                  compute_lm_ci = do_ci,
                                  Mboot = Mboot),
                 silent = FALSE)

    if (inherits(res_m, "try-error")) {
      cat(sprintf("  ERROR: replication %d failed, skipped.\n", m))
    } else {
      for (w in w_list) {
        if (!is.null(res_m[[as.character(w)]])) {
          results[[as.character(w)]][m, ] <- res_m[[as.character(w)]]
        }
      }
    }

    iter_dur <- as.numeric(difftime(Sys.time(), iter_start, units = "mins"))
    elapsed  <- as.numeric(difftime(Sys.time(), total_start, units = "mins"))
    avg_dur  <- elapsed / (m - start_idx + 1)
    remaining <- avg_dur * (M_total - m)
    cat(sprintf("  Time: %.1f min | Remaining: %.1f min (%.1f h)\n",
                iter_dur, remaining, remaining / 60))

    # SAVE CHECKPOINT
    if (m %% checkpoint_every == 0 || m == M_total) {
      saveRDS(list(results = results, next_idx = m + 1,
                   s = s, w_list = w_list, seed = seed,
                   timestamp = Sys.time()),
              file = checkpoint_file)
      cat(sprintf("  >>> Checkpoint saved (m=%d)\n", m))
    }
  }

  # Write the final CSVs
  for (w in w_list) {
    csv_path <- file.path(output_dir, sprintf("results_s%g_w%g.csv", s, w))
    write.csv(results[[as.character(w)]], csv_path, row.names = FALSE)
    cat(sprintf("CSV written: %s\n", csv_path))
  }

  total_time <- as.numeric(difftime(Sys.time(), total_start, units = "hours"))
  cat(sprintf("\n==== Completed (total %.1f h) ====\n", total_time))

  results
}

cat("run_sim_v6_core.R loaded.\n")
cat("Usage: run_scenario_v6(s = ..., w_list = c(...), output_dir = ...)\n")
