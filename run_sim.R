## ---- run_sim_v2.R ----

rm(list = ls()) # .rs.restartR()
Sys.setlocale(category = "LC_ALL", locale = "pt_PT.UTF-8")
pkgs <- c("survival","JMbayes2","nlme","lme4","riskRegression","prodlim","MASS","dplyr","splines","timeROC","haven")
need <- setdiff(pkgs, rownames(installed.packages()))
if(length(need)) install.packages(need, repos="https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only=TRUE))

source("survfitLMCR.R")

if (!dir.exists("results")) dir.create("results", recursive = TRUE)


# 1) Define the true parameters
true_params <- list(
  truebeta  = c("(Intercept)" = 10.63765,
                "year"        =  0.2629633,
                "drug"        = -0.09780603,
                "year:drug"   = -0.01708882),
  sigma     = 1.058143,  # residual SD
  Dmat      = matrix(c(0.6152592, 0.007199312, 0.007199312, 0.09915757),
                     nrow = 2, dimnames = list(c("(Intercept)","year"),
                                               c("(Intercept)","year"))),
  truegamma = c(transplanted = -0.4647676,
                dead         = 0.3482824),
  truealpha = c(transplanted = 0.1553865,
                dead         = 0.2506399)
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
n_intervals <- 1000 # Increased grid size from 200 to 1000 to avoid underestimation.
cuts    <- seq(0, T_max, length.out = n_intervals + 1)
delta   <- diff(cuts)
midpt   <- (cuts[-1] + cuts[-length(cuts)])/2

# Baseline hazard multipliers
f_tr <- 0.014  
f_de <- 0.0132


# spline-based baseline hazards
h_base_tr <- f_tr * exp(Bspline_basis(midpt) %*% gamma0)
h_base_de <- f_de * exp(Bspline_basis(midpt) %*% gamma0)

# 2a) Simulation function
simulate_joint_data_realistic <- function(n, params, CensorTime) {
  
  # --- Parameters ---
  q <- nrow(params$Dmat)
  b_mat <- MASS::mvrnorm(n, mu = rep(0, q), Sigma = params$Dmat)
  
  # Helper function: Event Time Simulation
  simulate_time_pc <- function(a_i, d_i, drug_i, gamma_k, alpha_k, h_base_cause, delta, midpt, cuts, CensorTime){
    U <- runif(1)
    H_target <- -log(U)
    lin_pred <- gamma_k * drug_i + alpha_k * (a_i + d_i * midpt)
    lin_pred <- pmin(lin_pred, 20) # Overflow protection. Increased limit from 5 to 20 to prevent JM overestimation.
    eff_haz  <- h_base_cause * exp(lin_pred)
    H_cum    <- cumsum(eff_haz * delta)
    j        <- which(H_cum >= H_target)[1]
    
    if (is.na(j)) return(CensorTime + 10) # If no event, push beyond the boundary
    
    H_prev <- if (j>1) H_cum[j-1] else 0
    t_prev <- cuts[j]
    t_prev + (H_target - H_prev) / eff_haz[j]
  }
  
  long_list  <- vector("list", n)
  surv_list  <- vector("list", n)
  truth_list <- vector("list", n)
  
  for (i in seq_len(n)) {
    
    # 1. Random Effects & Linear Predictors
    a_i <- params$truebeta["(Intercept)"] + b_mat[i,1] + params$covariates$drug[i] * params$truebeta["drug"]
    d_i <- params$truebeta["year"] + b_mat[i,2] + params$covariates$drug[i] * params$truebeta["year:drug"]
    drug_i <- params$covariates$drug[i]
    
    # 2. Competing Risks (Event Times)
    T1 <- simulate_time_pc(a_i, d_i, drug_i,
                           gamma_k = params$truegamma["transplanted"],
                           alpha_k = params$truealpha["transplanted"],
                           h_base_cause = h_base_tr,
                           delta, midpt, cuts, CensorTime)
    
    T2 <- simulate_time_pc(a_i, d_i, drug_i,
                           gamma_k = params$truegamma["dead"],
                           alpha_k = params$truealpha["dead"],
                           h_base_cause = h_base_de,
                           delta, midpt, cuts, CensorTime)
    
    T_true_event <- min(T1, T2)
    true_status  <- if (T1 < T2) 1 else 2
    
    # REALISTIC CENSORING (Exponential Drop-out) ---
    # Administrative censoring (15 years) + Random dropout (Exponential)
    Cens_admin  <- CensorTime
    Cens_random <- rexp(1, rate = 0.05) # Increasing the rate increases dropout
    Cens_i      <- min(Cens_admin, Cens_random)
    
    T_i     <- min(T_true_event, Cens_i)
    event_i <- if (Cens_i < T_true_event) 0 else true_status
    
    # VISIT TIMES (Jitter + Missing Visits) ---
    # Protocol: Baseline (0) followed by visits every year (or 0.8 years).
    protocol_visits <- seq(0, 15, by = 0.8) # E.g.: Every 0.8 years
    
    # Keep only times while subject is alive/under follow-up
    possible_visits <- protocol_visits[protocol_visits <= T_i]
    
    # If no visits remain (T_i is very small), at least include time 0
    if(length(possible_visits) == 0) possible_visits <- 0
    
    # A) Missing Visits: Probability of a visit occurring is 85%
    if (length(possible_visits) > 1) {
      kept_indices <- rbinom(length(possible_visits)-1, 1, 0.85) # Exclude time 0
      actual_visits <- c(0, possible_visits[-1][kept_indices == 1])
    } else {
      actual_visits <- possible_visits
    }
    
    # B) Jitter (Noise): Patients do not arrive exactly on schedule; introduce deviation.
    # sd = 0.1 (approx. 1 month deviation)
    t_ij <- actual_visits + rnorm(length(actual_visits), mean = 0, sd = 0.1)
    
    # Time cannot be negative or exceed T_i; baseline must be 0.
    t_ij[1] <- 0 
    t_ij <- sort(t_ij)
    t_ij <- t_ij[t_ij <= T_i]
    t_ij <- pmax(0, t_ij) # Set negative times to 0
    
    # --- Longitudinal Data Generation ---
    mu_ij <- a_i + d_i * t_ij
    Y_ij  <- mu_ij + rnorm(length(t_ij), 0, params$sigma)
    
    # Record
    long_list[[i]] <- data.frame(id=i, years=T_i, time=t_ij, Y=Y_ij, drug=drug_i)
    surv_list[[i]] <- data.frame(id=i, years=T_i, time=min(t_ij), Y=Y_ij[1], event=event_i, drug=drug_i)
    truth_list[[i]] <- data.frame(id=i, a_i=a_i, d_i=d_i)
  }
  
  long_df  <- do.call(rbind, long_list)
  surv_df  <- do.call(rbind, surv_list)
  truth_df <- do.call(rbind, truth_list)
  
  list(long=long_df, surv=surv_df, truth=truth_df)
}

# 2b) Monte-carlo simulation

simulate_once <- function(n, true_params, CensorTime, s = 2.5, horizon = 2, Mboot = 10){ 
  
  params <- true_params
  params$covariates <- list(drug = rbinom(n, 1, 0.5))
  
  sim <- simulate_joint_data_realistic(n, params, CensorTime)
  
  sim.surv  <- sim$surv
  sim.long  <- sim$long
  sim.truth <- sim$truth
  
  # table(sim.surv$event)
  # prop.table(table(sim.surv$event))
  # prop.table(table(pbc2.id$status))
  
  # Joint Model Setup
  simCR <- JMbayes2::crisk_setup(
    data        = sim.surv,
    statusVar   = "event",    # 0 = censored, 1 = transplanted, 2 = dead
    censLevel   = 0,
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
    data   = sim.long,
    control = lmeControl(opt = "optim") # Adding this is safe as optimization errors may occur with sparse data.
  )
  
  CR_forms <- list("Y" = ~ value(Y) : CR)
  
  jFit_CR <- JMbayes2::jm(
    Surv_object      = CoxFit_CR,
    Mixed_objects    = list(fm_Y),
    time_var         = "time",
    functional_forms = CR_forms,
    id_var           = "id",
    n_iter = 4000L,  
    n_burnin = 1000L,  
    n_thin = 3L,
    cores = 1 # Parallel processing might hide errors; set to 1 for testing
  )
  
  # summary(jFit_CR)
  
  
  
  
  
  
  # Joint Model Prediction
  t_star  <- s + horizon
  Ri      <- unique(sim.surv$id[ sim.surv$years > s ])
  
  # Expand dataframe (Allocate space for confidence intervals)
  jm_df <- data.frame(id = Ri, 
                      joint_trans = NA, trans_low = NA, trans_upp = NA,
                      joint_death = NA, death_low = NA, death_upp = NA)
  
  for(j in seq_along(Ri)){
    idj <- Ri[j]
    
    # Prepare patient data 
    ND_long  <- sim.long[sim.long$id == idj & sim.long$time < s, ]
    ND_event <- simCR[simCR$id == idj, ]
    ND_event$status2 <- 0; ND_event$years <- s
    ND <- list(newdataL = ND_long, newdataE = ND_event)
    
    # Get prediction
    P <- predict(jFit_CR,
                 newdata = ND,
                 return_newdata = TRUE,
                 times = t_star,
                 process = "event")
    
    
    # STRATA 2 = Transplant 
    jm_df$joint_trans[j] <- with(P, P$pred_CIF[`_strata`==2 & years == t_star])
    jm_df$trans_low[j]   <- with(P, P$low_CIF[`_strata`==2 & years == t_star])
    jm_df$trans_upp[j]   <- with(P, P$upp_CIF[`_strata`==2 & years == t_star])
    
    # STRATA 1 = Death
    jm_df$joint_death[j] <- with(P, P$pred_CIF[`_strata`==1 & years == t_star])
    jm_df$death_low[j]   <- with(P, P$low_CIF[`_strata`==1 & years == t_star])
    jm_df$death_upp[j]   <- with(P, P$upp_CIF[`_strata`==1 & years == t_star])
  }
  
  
  
  
  
  # TWO-STAGE LANDMARK APPROACH
  
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
    
    # Estimation of the linear mixed model
    lmeFit <- try(lme(fixed = Y ~ time*X, 
                      data = LMlong, 
                      random = ~ time | id, 
                      control = lmeControl(opt="optim", msMaxIter=1000, msMaxEval = 1000)), silent=TRUE)
    
    # If slope model fails, revert to simple model:
    if(inherits(lmeFit, "try-error")){
      lmeFit <- lme(fixed = Y ~ time*X, 
                    data = LMlong, 
                    random = ~ 1 | id,  # Random Intercept only
                    control = lmeControl(opt="optim"))
    }
    
    # BLUPs and parameters
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
    
    # Administrative censoring at the end of the prediction window
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
    
    ## Computation of the predicted individual cumulative incidences of events 
    # Subjects at risk at tLM
    Ri_pred <- LMsurv$id[which(LMsurv$years > tLM)]
    nLM_pred <- length(Ri_pred)
    # Computation of the predicted individual risk scores
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
    
    # Baseline hazards estimates
    bh1 <- basehaz(coxFit1, centered = FALSE)
    bh2 <- basehaz(coxFit2, centered = FALSE)
    
    # Reasonably quick function that converts cause-specific hazards to cumulative
    # incidence functions and extracts value at horizon
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
    
    # Individual cumulative incidences
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
  
  # TRUE landmark-conditional CIF (piecewise-constant exact within interval)
  true_cif_cond_at <- function(a_i, d_i, drug_i, s, t_star,
                               gamma_tr, alpha_tr, gamma_de, alpha_de,
                               h_base_tr, h_base_de, midpt, delta, cuts,
                               lin_cap = 5) {
    
    stopifnot(t_star >= s)
    
    # hazards on ALL intervals (defined at midpoints)
    lin1_all <- gamma_tr * drug_i + alpha_tr * (a_i + d_i * midpt)
    lin2_all <- gamma_de * drug_i + alpha_de * (a_i + d_i * midpt)
    lin1_all <- pmin(lin1_all, lin_cap)
    lin2_all <- pmin(lin2_all, lin_cap)
    
    h1_all   <- as.numeric(h_base_tr) * exp(lin1_all)
    h2_all   <- as.numeric(h_base_de) * exp(lin2_all)
    hsum_all <- h1_all + h2_all
    
    # helper: evaluate (S(t), CIF1(t), CIF2(t)) at any t in [0, max(cuts)]
    eval_at <- function(t) {
      if (t <= 0) return(list(S = 1, CIF1 = 0, CIF2 = 0))
      if (t >= max(cuts)) t <- max(cuts)
      
      # interval index j where cuts[j] < t <= cuts[j+1]
      j <- findInterval(t, cuts, rightmost.closed = TRUE)
      # findInterval returns length(cuts) when t==max(cuts); cap to last interval
      j <- min(j, length(delta))
      dt <- t - cuts[j]  # time within interval j
      
      # cumulative quantities up to start of interval j
      if (j == 1) {
        H_before    <- 0
        CIF1_before <- 0
        CIF2_before <- 0
        S_before    <- 1
      } else {
        H_inc    <- hsum_all * delta
        H_end    <- cumsum(H_inc)
        
        # exact CIF accumulation up to interval ends:
        # For constant hazards in interval j:
        # increment CIF1 over full interval = S_before * (h1/hsum) * (1 - exp(-hsum*delta))
        # similarly CIF2
        # We build those increments cumulatively.
        S_start <- c(1, exp(-H_end[-length(H_end)]))  # survival at interval starts
        
        inc1 <- S_start * (h1_all/hsum_all) * (1 - exp(-hsum_all * delta))
        inc2 <- S_start * (h2_all/hsum_all) * (1 - exp(-hsum_all * delta))
        
        CIF1_end <- cumsum(inc1)
        CIF2_end <- cumsum(inc2)
        
        H_before    <- H_end[j-1]
        CIF1_before <- CIF1_end[j-1]
        CIF2_before <- CIF2_end[j-1]
        S_before    <- exp(-H_before)
      }
      
      # now add partial interval contribution up to time t (dt within interval j)
      hz  <- hsum_all[j]
      if (hz <= 0) {
        # degenerate (shouldn't happen if hazards positive)
        H_t   <- H_before
        S_t   <- S_before
        CIF1_t <- CIF1_before
        CIF2_t <- CIF2_before
      } else {
        # exact partial increment with constant hazards:
        # ∫_0^dt S_before*exp(-hz x) * h1 dx = S_before * h1/hz * (1 - exp(-hz dt))
        fac <- (1 - exp(-hz * dt))
        CIF1_t <- CIF1_before + S_before * (h1_all[j]/hz) * fac
        CIF2_t <- CIF2_before + S_before * (h2_all[j]/hz) * fac
        H_t    <- H_before + hz * dt
        S_t    <- exp(-H_t)
      }
      
      list(S = S_t, CIF1 = CIF1_t, CIF2 = CIF2_t)
    }
    
    # evaluate at s and t_star
    at_s  <- eval_at(s)
    at_ts <- eval_at(t_star)
    
    # conditional (given event-free at s)
    S_s <- at_s$S
    if (S_s <= 0) {
      return(c(CIF1 = NA_real_, CIF2 = NA_real_))
    }
    
    CIF1_cond <- (at_ts$CIF1 - at_s$CIF1) / S_s
    CIF2_cond <- (at_ts$CIF2 - at_s$CIF2) / S_s
    
    # numerical safety
    CIF1_cond <- max(0, min(1, CIF1_cond))
    CIF2_cond <- max(0, min(1, CIF2_cond))
    if (CIF1_cond + CIF2_cond > 1) {
      # clip very small overshoots
      tot <- CIF1_cond + CIF2_cond
      CIF1_cond <- CIF1_cond / tot
      CIF2_cond <- CIF2_cond / tot
    }
    
    c(CIF1 = CIF1_cond, CIF2 = CIF2_cond)
  }
  
  
  # attach drug to truth
  truth_sub <- merge(
    sim.truth,
    sim.surv[, c("id", "drug")],
    by = "id",
    all.x = TRUE,
    sort = FALSE
  )
  truth_sub <- truth_sub[match(out$id, truth_sub$id), ]
  stopifnot(all(truth_sub$id == out$id))
  
  # parameters (explicit)
  gamma_tr <- params$truegamma["transplanted"]
  alpha_tr <- params$truealpha["transplanted"]
  gamma_de <- params$truegamma["dead"]
  alpha_de <- params$truealpha["dead"]
  
  # compute truth for each subject (matrix 2 x n)
  true_mat <- mapply(
    FUN = function(a, d, drug) {
      true_cif_cond_at(
        a_i = a, d_i = d, drug_i = drug,
        s = s, t_star = t_star,
        gamma_tr = gamma_tr, alpha_tr = alpha_tr,
        gamma_de = gamma_de, alpha_de = alpha_de,
        h_base_tr = h_base_tr, h_base_de = h_base_de,
        midpt = midpt, delta = delta, cuts = cuts,
        lin_cap = 5
      )
    },
    a    = truth_sub$a_i,
    d    = truth_sub$d_i,
    drug = truth_sub$drug
  )
  
  true_trans   <- as.numeric(true_mat["CIF1", ])
  true_death <- as.numeric(true_mat["CIF2", ])
  
  
  # Add TRUE values to the out dataframe
  out$TRUE_trans <- true_trans
  out$TRUE_death <- true_death
  
  
  
  # COVERAGE CALCULATION
  
  # 1. Transplant Coverage
  # Logic: If Lower Bound <= True Value <= Upper Bound then 1, else 0
  cov_trans_vec <- (true_trans >= jm_df$trans_low) & (true_trans <= jm_df$trans_upp)
  cov_trans_jm  <- mean(cov_trans_vec, na.rm = TRUE)
  
  # 2. Death Coverage
  cov_death_vec <- (true_death >= jm_df$death_low) & (true_death <= jm_df$death_upp)
  cov_death_jm  <- mean(cov_death_vec, na.rm = TRUE)
  
  
  
  
  
  
  
  stopifnot(all(is.finite(true_trans)), all(true_trans >= 0), all(true_trans <= 1))
  stopifnot(all(is.finite(true_death)), all(true_death >= 0), all(true_death <= 1))
  stopifnot(all(true_trans + true_death <= 1 + 1e-6))
  
  
  
  cat("\n--- DEBUG: ranges ---\n")
  print(range(out$joint_trans, na.rm=TRUE))
  print(range(out$joint_death, na.rm=TRUE))
  print(range(out$Rec_Mean, na.rm=TRUE))
  print(range(out$Death_Mean, na.rm=TRUE))
  
  cat("\n--- DEBUG: cor with truth ---\n")
  print(cor(out$joint_trans, true_trans, use="complete.obs"))
  print(cor(out$joint_death, true_death, use="complete.obs"))
  
  
  
  
  # Bias/RMSE
  bias_trans_jm    <- mean(out$joint_trans - true_trans)
  rmse_trans_jm    <- sqrt(mean((out$joint_trans - true_trans)^2))
  bias_death_jm  <- mean(out$joint_death - true_death)
  rmse_death_jm  <- sqrt(mean((out$joint_death - true_death)^2))
  
  bias_trans_lm    <- mean(out$Rec_Mean    - true_trans)
  rmse_trans_lm    <- sqrt(mean((out$Rec_Mean    - true_trans)^2))
  bias_death_lm  <- mean(out$Death_Mean  - true_death)
  rmse_death_lm  <- sqrt(mean((out$Death_Mean  - true_death)^2))
  
  # IPCW-corrected Brier & Time-dependent AUC
  # Evaluation data:
  evaldat <- sim.surv[match(out$id, sim.surv$id), c("years","event")]
  
  sc_trans_brier <- Score(
    object      = list(JM = out$joint_trans, LM = out$Rec_Mean),
    formula     = Hist(years, event) ~ 1,
    data        = evaldat,
    times       = t_star,
    cause       = 1,
    metrics     = "brier",
    cens.model = "km",
    conservative = TRUE
  )
  sc_de_brier <- Score(
    object      = list(JM = out$joint_death, LM = out$Death_Mean),
    formula     = Hist(years, event) ~ 1,
    data        = evaldat,
    times       = t_star,
    cause       = 2,
    metrics     = "brier",
    cens.model = "km",
    conservative = TRUE
  )
  brier_trans_jm   <- subset(sc_trans_brier$Brier$score,  model=="JM")$Brier[1]
  brier_trans_lm   <- subset(sc_trans_brier$Brier$score,  model=="LM")$Brier[1]
  brier_death_jm <- subset(sc_de_brier$Brier$score,   model=="JM")$Brier[1]
  brier_death_lm <- subset(sc_de_brier$Brier$score,   model=="LM")$Brier[1]
  
  sc_trans_auc <- Score(
    object      = list(JM = out$joint_trans, LM = out$Rec_Mean),
    formula     = Hist(years, event) ~ 1,
    data        = evaldat,
    times       = t_star,
    cause       = 1,
    metrics     = "auc",
    cens.model = "km",
    conservative = TRUE
  )
  sc_de_auc <- Score(
    object      = list(JM = out$joint_death, LM = out$Death_Mean),
    formula     = Hist(years, event) ~ 1,
    data        = evaldat,
    times       = t_star,
    cause       = 2,
    metrics     = "auc",
    cens.model = "km",
    conservative = TRUE
  )
  auc_trans_jm   <- subset(sc_trans_auc$AUC$score, model=="JM")$AUC[1]
  auc_trans_lm   <- subset(sc_trans_auc$AUC$score, model=="LM")$AUC[1]
  auc_death_jm <- subset(sc_de_auc$AUC$score,  model=="JM")$AUC[1]
  auc_death_lm <- subset(sc_de_auc$AUC$score,  model=="LM")$AUC[1]
  
  
  # Calculate Event Counts ---
  # Window of interest: (s, s + horizon]
  t_window_end <- s + horizon
  
  # Events only among those in the risk set (alive at time s)
  # (The condition years > s ensures this)
  n_event_trans_window <- sum(sim.surv$years > s & sim.surv$years <= t_window_end & sim.surv$event == 1)
  n_event_death_window <- sum(sim.surv$years > s & sim.surv$years <= t_window_end & sim.surv$event == 2)
  
  # ---------------------------------------------------
  
  # =========================================================================
  # TRANSFER DATA TO GLOBAL ENVIRONMENT
  # =========================================================================
  
  if(exists("jFit_CR")) jFit_CR <<- jFit_CR
  if(exists("P"))       P       <<- P
  if(exists("LMout_det")) LMout_det <<- LMout_det
  if(exists("out"))     out     <<- out
  if(exists("sim.long")) sim.long <<- sim.long 
  if(exists("sim.surv")) sim.surv <<- sim.surv
  if(exists("jm_df")) jm_df <<- jm_df
  
  # =========================================================================
  
  
  # Return: 21 metrics (bias/rmse JM & LM; IPCW-Brier; AUC)
  return(list(
    metrics = c(bias_trans_jm, rmse_trans_jm,
                bias_death_jm, rmse_death_jm,
                bias_trans_lm,  rmse_trans_lm,
                bias_death_lm,rmse_death_lm,
                brier_trans_jm, brier_trans_lm,
                brier_death_jm, brier_death_lm,
                auc_trans_jm, auc_trans_lm,
                auc_death_jm, auc_death_lm,
                cov_trans_jm,
                cov_death_jm,
                nrow(jm_df),            # Risk Set size
                n_event_trans_window,  # Number of transplants in window
                n_event_death_window), # Number of deaths in window
    out = out,
    truth = data.frame(id = out$id, true_trans = true_trans, true_death = true_death),
    evaldat = evaldat
  ))
}





# 3) Monte Carlo replications
M          <- 200
n          <- 500 # Sample size (500 or 1000)
CensorTime <- 15
params     <- true_params

results <- matrix(NA_real_, nrow = M, ncol = 21,
                  dimnames = list(NULL,
                                  c("bias_trans_jm", "rmse_trans_jm",
                                    "bias_death_jm", "rmse_death_jm",
                                    "bias_trans_lm", "rmse_trans_lm",
                                    "bias_death_lm", "rmse_death_lm",
                                    "brier_trans_jm", "brier_trans_lm",
                                    "brier_death_jm", "brier_death_lm",
                                    "auc_trans_jm", "auc_trans_lm",
                                    "auc_death_jm", "auc_death_lm",
                                    "cov_trans_jm", "cov_death_jm",
                                    "n_at_risk",
                                    "n_events_trans", "n_events_death")))


set.seed(2025)

# Get overall start time before loop
total_start_time <- Sys.time()

for(m in seq_len(M)){
  
  # 1. Record start time for this iteration
  iter_start_time <- Sys.time()
  
  
  cat("Simulation:", m, "/", M, "\n")
  res_m <- try(simulate_once(n, params, CensorTime,
                             s = 2.5, horizon = 2, Mboot = 10), 
               # s landmark time: (2.5, 5.5), horizon:(2, 4, 6)
               silent = TRUE)
  
  
  if(inherits(res_m, "try-error")){
    warning("Replication ", m, " failed; skipped.")
    # Calculate duration even if error occurs to track time loss
  } else {
    results[m, ] <- res_m$metrics
  }
  
  
  # 2. Record end time for this iteration and get difference
  iter_end_time <- Sys.time()
  iter_duration <- difftime(iter_end_time, iter_start_time, units = "mins")
  
  # 3. Estimate remaining time
  # Total time elapsed / completed iterations = Average iteration time
  avg_duration <- difftime(iter_end_time, total_start_time, units = "mins") / m
  remaining_iters <- M - m
  est_remaining_time <- avg_duration * remaining_iters
  
  # 4. Print Info to Screen
  cat(sprintf("[Done] Duration: %.2f min. | Est. Remaining: %.2f min. (%.1f hours)\n", 
              iter_duration, est_remaining_time, est_remaining_time/60))
  
}


results




# Inspect results
print(colMeans(results, na.rm = TRUE))


# Monte Carlo summary table
col_means <- colMeans(results, na.rm = TRUE)
col_sds   <- apply(results, 2, sd, na.rm = TRUE)

n_eff <- sum(complete.cases(results))
col_ses <- col_sds / sqrt(n_eff)

summary_tab <- data.frame(
  metric = colnames(results),
  mean   = as.numeric(col_means),
  sd     = as.numeric(col_sds),
  se     = as.numeric(col_ses),
  row.names = NULL
)
print(summary_tab, digits = 4)


# Save Results
if (!dir.exists("results")) dir.create("results", recursive = TRUE)
out_path <- file.path("results", "mc_summary_metrics_01_n500_s2.5_w2.csv")
write.csv(summary_tab, out_path, row.names = FALSE)
cat("\nSaved to:", normalizePath(out_path), "\n")

save.image(file = file.path("results", "simulation_workspace_01_n500_s2.5_w2.RData"))

write.csv(results, "results/results_01_n500_s2.5_w2_raw.csv", row.names = FALSE)
