## ---- run_sim.R ----
# Gerekli paketler
pkgs <- c("survival","JMbayes2","nlme","lme4","riskRegression","prodlim","MASS","dplyr","splines","timeROC","haven")
need <- setdiff(pkgs, rownames(installed.packages()))
if(length(need)) install.packages(need, repos="https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only=TRUE))

source("survfitLMCR.R")

if (!dir.exists("results")) dir.create("results", recursive = TRUE)

rm(list = ls()) # .rs.restartR()

Sys.setlocale(category = "LC_ALL", locale = "pt_PT.UTF-8")


# 1) Define the true parameters
true_params <- list(
  truebeta  = c("(Intercept)" = 10.64239343,
                "year"        =  0.24400503,
                "drug"        = -0.09458877,
                "year:drug"   = -0.01220294),
  sigma     = 1.058395,  # residual SD
  Dmat      = matrix(c(0.588910, -0.014538, -0.014538, 0.101550),
                     nrow = 2, dimnames = list(c("(Intercept)","year"),
                                               c("(Intercept)","year"))),
  truegamma = c(transplanted = 0.1553865,
                dead         = 0.2506399),
  alpha = 0.05 # longitudinal–survival association coefficient (0.05, 0.3, 0.5)
)

library(splines)

# spline-based baseline hazard parameters
T_max      <- 15
n_knots    <- 5
knots      <- seq(0, T_max, length = n_knots + 2)[-c(1, n_knots+2)]
Bspline_basis <- function(t) bs(t,
                                knots  = knots,
                                degree = 3,
                                intercept = TRUE)

# spline coefficients (example profile)
K       <- length(knots) + 3 + 1
gamma0 <- seq(-3, 0, length.out = K)

# piecewise grid
n_intervals <- 200
cuts    <- seq(0, T_max, length.out = n_intervals + 1)
delta   <- diff(cuts)
midpt   <- (cuts[-1] + cuts[-length(cuts)])/2

# Baseline hazard multipliers
f_tr <- 0.24
f_de <- 0.08

# spline-based baseline hazards
h_base_tr <- f_tr * exp(Bspline_basis(midpt) %*% gamma0)
h_base_de <- f_de * exp(Bspline_basis(midpt) %*% gamma0)

# 2) Simulation function
simulate_joint_data_exact <- function(n, params, times_empirical, CensorTime, alpha=0.5) {
  
  # subject-level covariate
  q     <- nrow(params$Dmat)
  b_mat <- MASS::mvrnorm(n, mu = rep(0, q), Sigma = params$Dmat)
  
  simulate_time_pc <- function(a_i, d_i, gamma_surv, h_base_cause, delta, midpt, cuts, CensorTime){
    U        <- runif(1)
    H_target <- -log(U)
    lin_pred <- gamma_surv * (a_i + alpha * d_i * midpt)
    lin_pred <- pmin(lin_pred, 5)
    eff_haz  <- h_base_cause * exp(lin_pred)
    H_cum    <- cumsum(eff_haz * delta)
    j        <- which(H_cum >= H_target)[1]
    if (is.na(j)) return(CensorTime)
    H_prev <- if (j>1) H_cum[j-1] else 0
    t_prev <- cuts[j]
    t_prev + (H_target - H_prev) / eff_haz[j]
  }
  
  long_list  <- vector("list", n)
  surv_list  <- vector("list", n)
  truth_list <- vector("list", n)
  
  for (i in seq_len(n)) {
    
    a_i <- params$truebeta["(Intercept)"] + b_mat[i,1] +
      params$covariates$drug[i] * params$truebeta["drug"]
    d_i <- params$truebeta["year"] + b_mat[i,2] +
      params$covariates$drug[i] * params$truebeta["year:drug"]
    
    # competing risks
    T1 <- simulate_time_pc(a_i, d_i,
                           params$truegamma["transplanted"],
                           h_base_tr, delta, midpt, cuts, CensorTime)
    
    T2 <- simulate_time_pc(a_i, d_i,
                           params$truegamma["dead"],
                           h_base_de, delta, midpt, cuts, CensorTime)
    
    T_event                       <- min(T1, T2)
    event_i                       <- if (T1 < T2) 1 else 2
    Cens_i                        <- runif(1, 0, CensorTime)
    T_i                           <- min(T_event, Cens_i)
    if (Cens_i < T_event) event_i <- 0
    
    # longitudinal
    m_i            <- sample(8:12,1)
    valid_times    <- times_empirical[times_empirical > 0 & times_empirical < T_i]
    if(length(valid_times) == 0){
      t_ij <- rep(0, m_i)
    } else {
      k    <- min(length(valid_times), m_i - 1)
      t_ij <- sort(c(0, sample(valid_times, size = k, replace = FALSE)))
    }
    k              <- min(length(valid_times), m_i - 1)
    t_ij           <- sort(c(0, sample(valid_times, size = k, replace = FALSE)))
    mu_ij          <- a_i + d_i * t_ij
    Y_ij           <- mu_ij + rnorm(length(t_ij), 0, params$sigma)
    
    long_list[[i]] <- data.frame(id=i, years=T_i, time=t_ij, Y=Y_ij,
                                 drug=params$covariates$drug[i])
    surv_list[[i]] <- data.frame(id=i, years=T_i, time=min(t_ij), Y = Y_ij[1],
                                 event=event_i, drug=params$covariates$drug[i])
    truth_list[[i]] <- data.frame(id=i, a_i=a_i, d_i=d_i)   }
  
  long_df  <- do.call(rbind, long_list);  rownames(long_df)  <- NULL
  surv_df  <- do.call(rbind, surv_list);  rownames(surv_df)  <- NULL
  truth_df <- do.call(rbind, truth_list); rownames(truth_df) <- NULL
  
  list(long=long_df, surv=surv_df, truth=truth_df)
}

# 2b) Monte-carlo simulation
simulate_once <- function(n, true_params, times_empirical, CensorTime, s = 5.5, horizon = 6, Mboot = 10){ 
  # s landmark time: (2.5, 5.5), horizon:(2, 4, 6)
  
  params <- true_params
  params$covariates <- list(drug = rbinom(n, 1, 0.5))
  sim <- simulate_joint_data_exact(n, params, times_empirical, CensorTime, alpha = params$alpha)
  
  sim.surv  <- sim$surv
  sim.long  <- sim$long
  sim.truth <- sim$truth
  
  # Joint Model
  simCR <- crisk_setup(
    data       = sim.surv,
    statusVar  = "event",    # 0 = censored, 1 = transplanted, 2 = dead
    censLevel  = 0,
    nameStrata = "CR"
  )
  
  CoxFit_CR <- coxph(
    Surv(years, status2) ~ drug * strata(CR),
    data = simCR,
    x    = TRUE,
    y    = TRUE
  )
  
  fm_Y <- lme(
    Y ~ time * drug,
    random = ~ time | id,
    data   = sim.long
  )
  
  CR_forms <- list("Y" = ~ value(Y) : CR)
  
  jFit_CR <- JMbayes2::jm(
    Surv_object      = CoxFit_CR,
    Mixed_objects    = list(fm_Y),
    time_var         = "time",
    functional_forms = CR_forms,
    id_var           = "id",
    n_iter = 3000L,
    n_burnin = 500L,
    n_thin = 3L
  )
  # summary(jFit_CR)
  
  # Joint Model Prediction
  t_star  <- s + horizon
  Ri      <- unique(sim.surv$id[ sim.surv$years > s ])
  
  jm_df <- data.frame(id=Ri, joint_trans=NA, joint_death=NA)
  for(j in seq_along(Ri)){
    idj <- Ri[j]
    ND_long  <- sim.long[sim.long$id==idj & sim.long$time < s, ]
    ND_event <- simCR[simCR$id==idj, ]
    ND_event$status2 <- 0; ND_event$years <- s
    ND <- list(newdataL = ND_long, newdataE = ND_event)
    
    P <- predict(jFit_CR,
                 newdata=ND,
                 return_newdata=TRUE,
                 times=t_star,
                 process="event")
    
    jm_df$joint_trans[j] <- with(P, P$pred_CIF[`_strata`==1 & years == t_star])
    jm_df$joint_death[j] <- with(P, P$pred_CIF[`_strata`==2 & years == t_star])
  }
  
  # TWO-STAGE LANDMARK APPROACH
  getwd()
  source("survfitLMCR.R")
  
  survfitLMCR <- function(tLM, thor, simulate = T, M = 1000, CI.levels = c(0.025, 0.975)) {
    tpred <- tLM + thor
    
    Ri     <- sim.surv$id[sim.surv$years > tLM]
    nLM <- length(Ri)
    LMlong <- subset(sim.long, time < tLM & id %in% Ri)
    LMlong <- LMlong %>%rename(X = drug)
    LMsurv <- subset(sim.surv, id %in% Ri)
    LMsurv <- LMsurv %>%rename(X = drug)
    LMsurv$Rec  <- as.integer(LMsurv$event == 1)
    LMsurv$Death <- as.integer(LMsurv$event == 2)
    
    lmeFit <- lme(fixed = Y ~ time*X,
                  data = LMlong,
                  random = ~ time| id,
                  control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))
    
    b <- ranef(lmeFit)
    sigma <- lmeFit$sigma
    D <- lapply(pdMatrix(lmeFit$modelStruct$reStruct), "*", sigma^2)[[1]]
    betas <- fixef(lmeFit)
    
    Data.s <- data.frame(id = unique(LMlong$id),
                         time = tLM,
                         X = LMlong[!duplicated(LMlong$id), "X"])
    Xlong.s <- model.matrix( ~ time*X, Data.s)
    Z.s <- model.matrix( ~ time, Data.s)
    LMsurv$level <- as.vector(c(Xlong.s %*% betas) + rowSums(Z.s * b))
    Xlong_deriv.s <- model.matrix( ~ X, Data.s)
    Zderiv.s <- model.matrix( ~ 1, Data.s)
    LMsurv$slope <- as.vector(c(Xlong_deriv.s %*% betas[c(2, 4)]) +
                                rowSums(Zderiv.s * b[ , 2, drop = FALSE]))
    
    LMsurv$tRecAC <- pmin(LMsurv$years, tpred)
    LMsurv$RecAC <- LMsurv$Rec
    LMsurv$RecAC[LMsurv$years > tpred] <- 0
    LMsurv$tDeathAC <- pmin(LMsurv$years, tpred)
    LMsurv$DeathAC <- LMsurv$Death
    LMsurv$DeathAC[LMsurv$years > tpred] <- 0
    
    coxFit1 <- coxph(Surv(tRecAC, RecAC) ~ X + level + slope, data = LMsurv,
                     x = TRUE, y = TRUE)
    coxFit2 <- coxph(Surv(tDeathAC, DeathAC) ~ X + level + slope, data = LMsurv,
                     x = TRUE, y = TRUE)
    
    Ri_pred <- LMsurv$id[which(LMsurv$years > tLM)]
    nLM_pred <- length(Ri_pred)
    LMsurv_pred <- LMsurv[LMsurv$id %in% Ri_pred, ]
    LMlong_pred <- LMlong[LMlong$time < tLM & LMlong$id %in% Ri_pred, ]
    
    Xlong_pred <- split(data.frame(model.matrix( ~ time*X, LMlong_pred)),LMlong_pred$id)
    Z_pred <- split(data.frame(model.matrix( ~ time, LMlong_pred)),
                    LMlong_pred$id)
    Y_pred <- split(LMlong_pred$Y, LMlong_pred$id)
    V_pred <- lapply(Z_pred, function(x){
      as.matrix(x) %*% D %*% t(as.matrix(x)) + diag(sigma^2, nrow(as.matrix(x)))})
    
    b_pred <- lapply(seq_len(length(Xlong_pred)),
                     function(i, x, y, z, v) {
                       D %*% t(as.matrix(z[[i]])) %*%
                         solve(v[[i]]) %*% (y[[i]] - as.matrix(x[[i]]) %*% fixef(lmeFit))},
                     x = Xlong_pred, y = Y_pred, z = Z_pred, v = V_pred)
    b_pred <- matrix(unlist(b_pred), ncol = 2, , byrow = T)
    Data_pred.s <- data.frame(id = unique(LMlong_pred$id),
                              time = tLM,
                              X = LMlong_pred[!duplicated(LMlong_pred$id), "X"])
    Xlong_pred.s <- model.matrix( ~ time*X, Data_pred.s)
    Z_pred.s <- model.matrix( ~ time, Data_pred.s)
    LMsurv_pred$level <- as.vector(c(Xlong_pred.s %*% betas) + rowSums(Z_pred.s * b_pred))
    Xlong_deriv_pred.s <- model.matrix( ~ X, Data_pred.s)
    Zderiv_pred.s <- model.matrix( ~ 1, Data_pred.s)
    LMsurv_pred$slope <- as.vector(c(Xlong_deriv_pred.s %*% betas[c(2, 4)]) +
                                     rowSums(Zderiv_pred.s * b_pred[ , 2, drop = FALSE]))
    Xsurv_pred <- model.matrix( ~ 0 + X + level + slope, LMsurv_pred)
    HR1 <- as.numeric(exp(Xsurv_pred %*% coxFit1$coef))
    HR2 <- as.numeric(exp(Xsurv_pred %*% coxFit2$coef))
    
    bh1 <- basehaz(coxFit1, centered = FALSE)
    bh2 <- basehaz(coxFit2, centered = FALSE)
    
    toci <- function(bh1, bh2, HR1, HR2, tpred)
    {
      h1 <- bh1
      names(h1)[1] <- "Haz"
      h1$Haz <- h1$Haz * HR1
      h1$cause <- 1
      h2 <- bh2
      names(h2)[1] <- "Haz"
      h2$Haz <- h2$Haz * HR2
      h2$cause <- 2
      Haz <- rbind(h1, h2)
      CI <- CumInc(Haz)
      idx <- sum(CI$time <= tpred)
      return(CI[idx, ])
    }
    
    ci <- matrix(NA, nLM_pred, 2)
    for (i in 1:nLM_pred) {
      ci[i, ] <- as.numeric(toci(bh1, bh2, HR1[i], HR2[i], tpred))[2:3]
      ci
    }
    
    res <- replicate(2, matrix(NA, nrow = nLM_pred, ncol = 3), simplify = FALSE)
    for(k in 1:2) {
      for (i in seq_len(nLM_pred)) {
        res[[k]][i, ] <- c(tLM, thor, ci[i, k])
      }
      colnames(res[[k]]) <- c("tLM", "thor", "Value")
      rownames(res[[k]]) <- sapply(Ri_pred, function(x) paste0("ID", x))
    }
    names(res) <- c("Rec", "Death")
    
    res.MC <- list(Rec = NULL, Death = NULL)
    names(res.MC) <- c("Rec", "Death")
    
    result <- list(res = res, res.MC = res.MC, simulate = FALSE, M = 0)
    if (exists(".Random.seed", envir = globalenv())) {
      rm(list = ".Random.seed", envir = globalenv())
    }
    class(result) <- "survfitLMCR"
    return(result)
  }
  
  # Landmark Model
  LMout_det <- survfitLMCR(tLM = s, thor = horizon, simulate = FALSE)
  Rec_df   <- as.data.frame(LMout_det$res$Rec);   Rec_df   <- Rec_df["Value"];   names(Rec_df)   <- "Rec_Mean"
  Death_df <- as.data.frame(LMout_det$res$Death); Death_df <- Death_df["Value"]; names(Death_df) <- "Death_Mean"
  lm_df <- merge(Rec_df, Death_df, by="row.names")
  names(lm_df)[1] <- "id"
  lm_df$id <- as.integer(gsub("ID *", "", lm_df$id))
  
  out <- merge(jm_df, lm_df, by="id")
  
  # TRUE CIF
  true_cif_at <- function(a_i, d_i, t_star){
    idx <- which(midpt <= t_star)
    lin1 <- params$truegamma["transplanted"] * (a_i + params$alpha * d_i * midpt[idx])
    lin2 <- params$truegamma["dead"]         * (a_i + params$alpha * d_i * midpt[idx])
    h1   <- h_base_tr[idx] * exp(pmin(lin1, 5))
    h2   <- h_base_de[idx] * exp(pmin(lin2, 5))
    hsum <- h1 + h2
    Hcum <- cumsum(hsum * delta[idx])
    S    <- exp(-Hcum)
    S_lag <- c(1, S[-length(S)])
    CIF1 <- sum(h1 * S_lag * delta[idx])
    CIF2 <- sum(h2 * S_lag * delta[idx])
    c(CIF1 = CIF1, CIF2 = CIF2)
  }
  truth_sub <- sim.truth[match(out$id, sim.truth$id), ]
  stopifnot(all(truth_sub$id == out$id))
  true_rec   <- mapply(function(a,d) true_cif_at(a,d,t_star)[["CIF1"]], truth_sub$a_i, truth_sub$d_i)
  true_death <- mapply(function(a,d) true_cif_at(a,d,t_star)[["CIF2"]], truth_sub$a_i, truth_sub$d_i)
  
  # Bias/RMSE
  bias_rec_jm    <- mean(out$joint_trans - true_rec)
  rmse_rec_jm    <- sqrt(mean((out$joint_trans - true_rec)^2))
  bias_death_jm  <- mean(out$joint_death - true_death)
  rmse_death_jm  <- sqrt(mean((out$joint_death - true_death)^2))
  
  bias_rec_lm    <- mean(out$Rec_Mean    - true_rec)
  rmse_rec_lm    <- sqrt(mean((out$Rec_Mean    - true_rec)^2))
  bias_death_lm  <- mean(out$Death_Mean  - true_death)
  rmse_death_lm  <- sqrt(mean((out$Death_Mean  - true_death)^2))
  
  # IPCW-corrected Brier & Time-dependent AUC
  # Evaluation data:
  evaldat <- sim.surv[match(out$id, sim.surv$id), c("years","event")]
  
  sc_rec_brier <- Score(
    object     = list(JM = out$joint_trans, LM = out$Rec_Mean),
    formula    = Hist(years, event) ~ 1,
    data       = evaldat,
    times      = t_star,
    cause      = 1,
    metrics    = "brier",
    cens.model = "km",
    conservative = TRUE
  )
  sc_de_brier <- Score(
    object     = list(JM = out$joint_death, LM = out$Death_Mean),
    formula    = Hist(years, event) ~ 1,
    data       = evaldat,
    times      = t_star,
    cause      = 2,
    metrics    = "brier",
    cens.model = "km",
    conservative = TRUE
  )
  brier_rec_jm   <- subset(sc_rec_brier$Brier$score,  model=="JM")$Brier[1]
  brier_rec_lm   <- subset(sc_rec_brier$Brier$score,  model=="LM")$Brier[1]
  brier_death_jm <- subset(sc_de_brier$Brier$score,   model=="JM")$Brier[1]
  brier_death_lm <- subset(sc_de_brier$Brier$score,   model=="LM")$Brier[1]
  
  sc_rec_auc <- Score(
    object     = list(JM = out$joint_trans, LM = out$Rec_Mean),
    formula    = Hist(years, event) ~ 1,
    data       = evaldat,
    times      = t_star,
    cause      = 1,
    metrics    = "auc",
    cens.model = "km",
    conservative = TRUE
  )
  sc_de_auc <- Score(
    object     = list(JM = out$joint_death, LM = out$Death_Mean),
    formula    = Hist(years, event) ~ 1,
    data       = evaldat,
    times      = t_star,
    cause      = 2,
    metrics    = "auc",
    cens.model = "km",
    conservative = TRUE
  )
  auc_rec_jm   <- subset(sc_rec_auc$AUC$score, model=="JM")$AUC[1]
  auc_rec_lm   <- subset(sc_rec_auc$AUC$score, model=="LM")$AUC[1]
  auc_death_jm <- subset(sc_de_auc$AUC$score,  model=="JM")$AUC[1]
  auc_death_lm <- subset(sc_de_auc$AUC$score,  model=="LM")$AUC[1]
  
  # Return: 16 metrics (bias/rmse JM & LM; IPCW-Brier; AUC)
  return(c(bias_rec_jm, rmse_rec_jm,
           bias_death_jm, rmse_death_jm,
           bias_rec_lm,  rmse_rec_lm,
           bias_death_lm,rmse_death_lm,
           brier_rec_jm, brier_rec_lm,
           brier_death_jm, brier_death_lm,
           auc_rec_jm, auc_rec_lm,
           auc_death_jm, auc_death_lm))
}

# 3) Monte Carlo replications
M          <- 200
n          <- 200 # Sample size (200 or 500)
times_empirical  <- sort(unique(pbc2$year))
CensorTime <- 15
params     <- true_params

results <- matrix(NA_real_, nrow = M, ncol = 16,
                  dimnames = list(NULL,
                                  c("bias_rec_jm","rmse_rec_jm",
                                    "bias_death_jm","rmse_death_jm",
                                    "bias_rec_lm","rmse_rec_lm",
                                    "bias_death_lm","rmse_death_lm",
                                    "brier_rec_jm","brier_rec_lm",
                                    "brier_death_jm","brier_death_lm",
                                    "auc_rec_jm","auc_rec_lm",
                                    "auc_death_jm","auc_death_lm")))

set.seed(2025)
for(m in seq_len(M)){
  cat("Simulation:", m, "/", M, "\n")
  res_m <- try(simulate_once(n, params, times_empirical, CensorTime,
                             s = 5.5, horizon = 6, Mboot = 10),  # s landmark time: (2.5, 5.5), horizon:(2, 4, 6)
               silent = TRUE)
  
  if(inherits(res_m, "try-error")){
    warning("Replication ", m, " failed; skipped.")
    next
  }
  results[m, ] <- res_m
}

results

# Monte Carlo summary table (meaningful when M > 1)
col_means <- colMeans(results, na.rm = TRUE)
col_sds   <- apply(results, 2, sd, na.rm = TRUE)
col_ses   <- col_sds / sqrt(nrow(results))

summary_tab <- data.frame(
  metric = colnames(results),
  mean   = as.numeric(col_means),
  sd     = as.numeric(col_sds),
  se     = as.numeric(col_ses),
  row.names = NULL
)
print(summary_tab, digits = 4)

# Optional: Save as CSV
# write.csv(summary_tab, "mc_summary_metrics.csv", row.names = FALSE)

