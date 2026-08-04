nwn <- function(x) suppressWarnings(as.numeric(x))

checkInputMatrix <- function(m, def, defDuration = NA) {
  m[,'dose'] <- nwn(m[,'dose'])
  m[!is.na(m[,'dose']) & m[,'dose'] < 0, 'dose'] <- NA
  if(!is.na(defDuration)) {
    if(!'duration' %in% names(m)) {
      m[,'duration'] <- defDuration
    } else {
      m[,'duration'] <- nwn(m[,'duration'])
      m[is.na(m[,'duration']), 'duration'] <- defDuration
      m[m[,'duration'] < 0, 'duration'] <- defDuration
    }
  }
  m[,'time'] <- nwn(m[,'dt']) 
  m <- m[complete.cases(m),]
  if(nrow(m) == 0) {
    m <- def
  }
  m
}

##########################
###### MY MODEL ##########
##########################
model_mtx <- function(usrbsa, usrSCR, usrgender) { #define model parameters 
  THETA1 <- 9.41216979510633  # CL pop
  THETA2 <- 30.666641775106   # V pop
  THETA3 <- 0.866697110621852 # Q2 pop
  THETA4 <- 5.61864647680956 # V2 pop
  THETA5 <- 0.13515317483686 # Q3 pop
  THETA6 <- 10.6049699903601 #V3 pop
  beta_scr_cl = -0.556395940439738 ; beta_bsa_cl = 0.609402364492555 ; beta_gender = 0.013242740163542
  
  sigprop <- 0.2508599^2
  sigadd = 0
  Omega = diag(c(0.242508391327, 0.289352674694, 0.38419083468, #cl v1 q2
                 0.300674583119, 0.466055297477, 0.343596438913)) #v2 q3 v3
  Sigma = diag(c(sigprop,sigadd))
  list(ncpt = 3, Omega = Omega, Sigma = Sigma,
       THETA1 = THETA1, THETA2 = THETA2, THETA3 = THETA3, THETA4 = THETA4, THETA5 = THETA5, THETA6 = THETA6,
       beta_scr_cl = beta_scr_cl, beta_bsa_cl = beta_bsa_cl, beta_gender = beta_gender,
       bsa = usrbsa, SCR_mgdl = usrSCR, pt_gender = usrgender,beta_scr_cl = beta_scr_cl, beta_bsa_cl = beta_bsa_cl, beta_gender = beta_gender )
}

mymodel <- function(time, state, parms) {
  with(as.list(c(state, parms)), {
    
    ## --- Interpolate covariates at current time ---
    bsa_t <- cov_bsa[!is.na(cov_bsa)][1]
    #bsa_t     <- approx(covtimes, cov_bsa,    xout = time, rule = 2)$y
    #SCR_mmol <- unique(cov_scr[!is.na(cov_scr)])
    SCR_fun <- approxfun(covtimes, cov_scr, method = "constant", rule = 2)
    SCR_mmol <- SCR_fun(time)    
    pt_gender <- ifelse(cov_gender[1] == "Male", 0 , 1)
    beta_bsa_cl = parms[["beta_bsa_cl"]]
    beta_scr_cl = parms[["beta_scr_cl"]]
    beta_gender = parms[["beta_gender"]]
    
    
    ## --- Structural parameters ---
    # You can include covariate effects here later (e.g., CL ∝ BSA^0.75)
    CL <- THETA1 * ((bsa_t/1.97)^beta_bsa_cl)*((SCR_mmol/68.08)^(beta_scr_cl))*exp(beta_gender*pt_gender) * exp(ETA1)
    V1 <- THETA2 * (bsa_t/1.97) * exp(ETA2)
    Q2 <- THETA3 * (bsa_t/1.97) * exp(ETA3)
    V2 <- THETA4 * (bsa_t/1.97) * exp(ETA4)
    Q3 <- THETA5 * (bsa_t/1.97) * exp(ETA5)
    V3 <- THETA6 * (bsa_t/1.97) * exp(ETA6)
    
    ## --- Microconstants ---
    k10 <- CL / V1
    k12 <- Q2 / V1
    k21 <- Q2 / V2
    k13 <- Q3 / V1
    k31 <- Q3 / V3
    
    ## --- Infusion rate ---
    rate <- 0
    for (i in seq_along(dose_times)) {
      if (time >= dose_times[i] && time <= (dose_times[i] + dose_durs[i])) {
        rate <- rate + dose_amts[i] / dose_durs[i]
      }
    }
    
    ## --- Differential equations ---
    dA1dt <- rate - (k10 + k12 + k13) * A1 + k21 * A2 + k31 * A3
    dA2dt <- k12 * A1 - k21 * A2
    dA3dt <- k13 * A1 - k31 * A3
    
    list(c(dA1dt, dA2dt, dA3dt)) 
  })
}

simulate_mymodel <- function(dat_subj, THETA1, THETA2, THETA3, THETA4, THETA5, THETA6, 
                             ETA1 = 0, ETA2 = 0, ETA3 = 0, ETA4 = 0, ETA5 = 0, ETA6 = 0, sigma = 0.01,
                             beta_bsa_cl, beta_scr_cl, beta_gender) {
  # ensure sorted
  dat_subj <- dat_subj %>% arrange(time)
  
  # extract doses
  dose_amts  <- dat_subj$amt[dat_subj$mdv == 1 & !is.na(dat_subj$amt)]
  dose_times <- dat_subj$time[dat_subj$mdv == 1 & !is.na(dat_subj$amt)]
  dose_durs  <- dat_subj$dur[dat_subj$mdv == 1 & !is.na(dat_subj$amt)]
  if (length(dose_durs) == 0) dose_durs <- rep(0, length(dose_amts))
  
  # pack parameters
  parms <- list(
    THETA1 = THETA1, THETA2 = THETA2, THETA3 = THETA3,
    THETA4 = THETA4, THETA5 = THETA5, THETA6 = THETA6,
    ETA1 = ETA1, ETA2 = ETA2, ETA3 = ETA3,
    ETA4 = ETA4, ETA5 = ETA5, ETA6 = ETA6,
    covtimes   = dat_subj$time,
    cov_bsa    = dat_subj$bsa,
    cov_scr    = dat_subj$SCR_mmol,
    cov_gender = dat_subj$pt_gender,
    dose_times = dose_times,
    dose_amts  = dose_amts,
    dose_durs  = dose_durs,
    beta_bsa_cl = beta_bsa_cl,
    beta_scr_cl = beta_scr_cl,
    beta_gender = beta_gender
  )
  
  # initial state
  state <- c(A1 = 0, A2 = 0, A3 = 0)
  print(max(dat_subj$time))
  times <- sort(unique(c(seq(0, max(dat_subj$time), by = 0.1), dat_subj$time)))
  
  # integrate
  #figure out tolerance
  out <- ode(y = state, times = times, func = mymodel, parms = parms, hmax = Inf)
  out <- as.data.frame(out)
  out$subject_id <- dat_subj$subject_id[1]
  
  SCR_fun <- approxfun(parms$covtimes, parms$cov_scr, method = "constant", rule = 2)
  SCR_mmol <- SCR_fun(times)    
  
  out$SCR_mmol = SCR_mmol
  
  # compute concentrations
  out <- out %>%
    rowwise() %>%
    mutate(
      bsa = dat_subj$bsa[1], #approx(parms$covtimes, parms$cov_bsa, xout = time, rule = 2)$y,
      pt_gender =  ifelse(dat_subj$pt_gender[1] == "Male", 0 , 1),
      CL = THETA1 * ((bsa/1.97)^beta_bsa_cl)*((SCR_mmol/68.08)^(beta_scr_cl))*exp(beta_gender*pt_gender)*exp(ETA1),
      V1 = THETA2 * (bsa/1.97)* exp(ETA2),
      Q2 = THETA3 * (bsa/1.97)* exp(ETA3),
      V2 = THETA4 * (bsa/1.97)* exp(ETA4),
      Q3 = THETA5 * (bsa/1.97)* exp(ETA5),
      V3 = THETA6 * (bsa/1.97)* exp(ETA6),
      Conc = A1 / V1
    ) %>%
    ungroup()
  
  out <- subset(out, time %in% dat_subj$time)
  out$evid <- dat_subj$evid
  
  return(out)
}

mapbayes_mymodel = function(eta,dat_subj, y, yt, Omega=mp$Omega, Sigma =mp$Sigma,
                            THETA1 = mp$THETA1, THETA2 = mp$THETA2, THETA3 = mp$THETA3, THETA4 = mp$THETA4, THETA5 = mp$THETA5, THETA6 = mp$THETA6,
                            beta_scr_cl = mp$beta_scr_cl, beta_bsa_cl = mp$beta_bsa_cl, beta_gender= mp$beta_gender){
  #suppressMessages(attach(mp))
  sigprop <- as.numeric(Sigma[1,1])
  sigadd <- as.numeric(Sigma[2,2])
  eta =  eta %>%  as.list
  # names(eta) <- names(init)
  eta_m <- eta %>% unlist %>% matrix(nrow=1)
  cat("mymodel eta:", paste(names(eta), round(unlist(eta), 4), sep="=", collapse=", "), "\n")
  
  looppat <- simulate_mymodel(dat_subj, THETA1, THETA2, THETA3, THETA4, THETA5, THETA6, 
                              ETA1 = eta$ETAcl, ETA2 = eta$ETAv1, 
                              ETA3=eta$ETAq2, ETA4=eta$ETAv2, 
                              ETA5 = eta$ETAq3, ETA6 = eta$ETAv3, beta_scr_cl = beta_scr_cl, beta_bsa_cl = beta_bsa_cl, beta_gender= beta_gender)
  
  newobs = looppat$Conc[which(looppat$time%in%yt)]
  
  # http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3339294/
  sig2j <- (newobs*sqrt(sigprop) + sqrt(sigadd))^2    #<------------------Check against pubmed
  sqwres <- log(sig2j) + (1/sig2j)*(y-newobs)^2
  
  nOn <- diag(eta_m %*% solve(Omega) %*% t(eta_m))
  return(sum(sqwres) + nOn)
}

pkprof_est_mymodel <- function(time_and_covs, mp) { 
  ncpt <- mp$ncpt
  
  y <- time_and_covs$dv[which(!is.na(time_and_covs$dv))[1]]
  yt <- round(time_and_covs$time[which(!is.na(time_and_covs$dv))[1]],1)
  
  init = eta <- c(ETAcl=.05, ETAv1=.05, ETAq2=.05, ETAv2=.05, ETAq3=.05, ETAv3=.05)
  #  print(time_and_covs$SCR_mmol)
  dat_subj = subset(time_and_covs, !is.na(dv) | !is.na(amt))
  print(dat_subj)
  fit <- newuoa(par=init,
                fn = mapbayes_mymodel,
                dat_subj = dat_subj,
                y=y,
                yt=yt,
                Omega=mp$Omega,
                Sigma =mp$Sigma,
                THETA1 = mp$THETA1, THETA2 = mp$THETA2, THETA3 = mp$THETA3, THETA4 = mp$THETA4, THETA5 = mp$THETA5, THETA6 = mp$THETA6,
                beta_scr_cl = mp$beta_scr_cl, beta_bsa_cl = mp$beta_bsa_cl, beta_gender= mp$beta_gender)
  
  return(fit$par)
}

##########################
######### HUI et al. #####
## ALL 2-compartment model
## Hui et al., J Clin Pharmacol 2019, 59(4):566-577
## Table 3: Final ALL model parameter estimates
## Covariates normalised by their sample medians:
##   BSA median = 0.735 m2, CLCR median = 192 mL/min/1.73m2, AGE median = 63.5 months
## CLCR computed internally via paediatric Schwartz formula:
##   CLCR (mL/min/1.73m2) = (0.413 * height_cm) / SCR_mgdl
##########################
model_hui <- function(usrbsa, usrSCR_mgdl, usr_height_cm, usr_age_months) {
  # --- Fixed effects (Table 3) ---
  TVCL  <- 7.73   # L/h
  TVV1  <- 19.0   # L
  TVQ   <- 0.283  # L/h
  TVV2  <- 6.63   # L
  
  # --- Covariate effect exponents ---
  theta_bsa_cl   <- 0.721   # BSA effect on CL
  theta_clcr_cl  <- 0.256   # CLCR effect on CL
  theta_bsa_v1   <- 0.985   # BSA effect on V1
  theta_age_q    <- 0.278   # AGE effect on Q
  
  # --- Covariate medians (Table 1, ALL cohort) ---
  med_bsa  <- 0.735   # m2
  med_clcr <- 192     # mL/min/1.73m2
  med_age  <- 63.5    # months
  
  # --- Compute CLCR via Schwartz formula ---
  # CLCR = 0.413 * height_cm / SCR_mgdl  (mL/min/1.73m2)
  clcr <- (0.413 * usr_height_cm) / usrSCR_mgdl
  
  # --- Omega (IIV as variances; Table 3 reports CV%) ---
  # ωCL = 14.3% CV  -> variance = log(1 + 0.143^2) ≈ 0.143^2 = 0.02045
  # ωIOV on CL = 14.9% -> same scale, treated as additional CL variability
  # ωV1 fixed at 0
  # ωQ  fixed at 0
  # ωV2 = 34.6% CV  -> variance ≈ 0.346^2 = 0.11972
  omega_cl  <- 0.143^2
  omega_v2  <- 0.346^2
  
  # 3 ETAs: ETAcl (IIV on CL), ETAiov (IOV on CL), ETAv2 (IIV on V2)
  Omega <- diag(c(omega_cl, omega_v2))
  
  # --- Residual (proportional, sigma = 30.2% CV) ---
  sigprop <- 0.302^2
  sigadd  <- 0
  Sigma   <- diag(c(sigprop, sigadd))
  
  list(
    ncpt          = 2,
    Omega         = Omega,
    Sigma         = Sigma,
    TVCL          = TVCL,
    TVV1          = TVV1,
    TVQ           = TVQ,
    TVV2          = TVV2,
    theta_bsa_cl  = theta_bsa_cl,
    theta_clcr_cl = theta_clcr_cl,
    theta_bsa_v1  = theta_bsa_v1,
    theta_age_q   = theta_age_q,
    med_bsa       = med_bsa,
    med_clcr      = med_clcr,
    med_age       = med_age,
    bsa           = usrbsa,
    SCR_mgdl      = usrSCR_mgdl,
    height_cm     = usr_height_cm,
    age_months    = usr_age_months,
    clcr          = clcr
  )
}

hui_ode <- function(time, state, parms) {
  with(as.list(c(state, parms)), {
    
    bsa_t <- cov_bsa[!is.na(cov_bsa)][1]
    
    SCR_fun  <- approxfun(covtimes, cov_scr_mgdl, method = "constant", rule = 2)
    scr_t    <- SCR_fun(time)
    clcr_t   <- (0.413 * height_cm) / scr_t
    
    CL <- TVCL * (bsa_t / med_bsa)^theta_bsa_cl * (clcr_t / med_clcr)^theta_clcr_cl * exp(ETAcl)
    V1 <- TVV1 * (bsa_t / med_bsa)^theta_bsa_v1
    Q  <- TVQ  * (age_months / med_age)^theta_age_q
    V2 <- TVV2 * exp(ETAv2)
    
    k10 <- CL / V1
    k12 <- Q  / V1
    k21 <- Q  / V2
    
    rate <- 0
    for (i in seq_along(dose_times)) {
      if (time >= dose_times[i] && time <= (dose_times[i] + dose_durs[i])) {
        rate <- rate + dose_amts[i] / dose_durs[i]
      }
    }
    
    dA1dt <- rate - (k10 + k12) * A1 + k21 * A2
    dA2dt <- k12 * A1 - k21 * A2
    
    list(c(dA1dt, dA2dt))
  })
}

simulate_hui <- function(dat_subj,
                         TVCL, TVV1, TVQ, TVV2,
                         theta_bsa_cl, theta_clcr_cl, theta_bsa_v1, theta_age_q,
                         med_bsa, med_clcr, med_age,
                         height_cm, age_months,
                         ETAcl = 0, ETAv2 = 0) {
  
  dat_subj <- dat_subj %>% arrange(time)
  
  dose_amts  <- dat_subj$amt[dat_subj$mdv == 1 & !is.na(dat_subj$amt)]
  dose_times <- dat_subj$time[dat_subj$mdv == 1 & !is.na(dat_subj$amt)]
  dose_durs  <- dat_subj$dur[dat_subj$mdv == 1 & !is.na(dat_subj$amt)]
  if (length(dose_durs) == 0) dose_durs <- rep(0, length(dose_amts))
  
  parms <- list(
    TVCL = TVCL, TVV1 = TVV1, TVQ = TVQ, TVV2 = TVV2,
    theta_bsa_cl  = theta_bsa_cl,
    theta_clcr_cl = theta_clcr_cl,
    theta_bsa_v1  = theta_bsa_v1,
    theta_age_q   = theta_age_q,
    med_bsa    = med_bsa,
    med_clcr   = med_clcr,
    med_age    = med_age,
    ETAcl      = ETAcl,
    ETAv2      = ETAv2,
    height_cm  = height_cm,
    age_months = age_months,
    covtimes      = dat_subj$time,
    cov_bsa       = dat_subj$bsa,
    cov_scr_mgdl  = dat_subj$SCR_mgdl,
    dose_times = dose_times,
    dose_amts  = dose_amts,
    dose_durs  = dose_durs
  )
  
  state <- c(A1 = 0, A2 = 0)
  times <- sort(unique(c(seq(0, max(dat_subj$time), by = 0.1), dat_subj$time)))
  
  out <- ode(y = state, times = times, func = hui_ode, parms = parms, hmax = Inf)
  out <- as.data.frame(out)
  out$subject_id <- dat_subj$subject_id[1]
  
  SCR_fun <- approxfun(parms$covtimes, parms$cov_scr_mgdl, method = "constant", rule = 2)
  scr_mgdl_vec <- SCR_fun(times)
  out$SCR_mmol <- scr_mgdl_vec * 88.4  # keep column name consistent with rest of app
  
  out <- out %>%
    rowwise() %>%
    mutate(
      bsa       = dat_subj$bsa[1],
      pt_gender = dat_subj$pt_gender[1],
      clcr_t    = (0.413 * height_cm) / SCR_fun(time),
      CL  = TVCL * (bsa / med_bsa)^theta_bsa_cl * (clcr_t / med_clcr)^theta_clcr_cl * exp(ETAcl),
      V1  = TVV1 * (bsa / med_bsa)^theta_bsa_v1,
      Q   = TVQ  * (age_months / med_age)^theta_age_q,
      V2  = TVV2 * exp(ETAv2),
      Conc = A1 / V1
    ) %>%
    ungroup()
  
  # Tolerance-based time match to avoid floating point mismatches
  tol <- 1e-6
  keep_idx <- sapply(dat_subj$time, function(t) {
    i <- which(abs(out$time - t) < tol)
    if (length(i) == 0) NA_integer_ else i[1]
  })
  valid <- !is.na(keep_idx)
  out <- out[keep_idx[valid], ]
  out$evid <- dat_subj$evid[valid]
  
  
  return(out)
}

mapbayes_hui <- function(eta, dat_subj, y, yt,
                         Omega = mp_hui$Omega, Sigma = mp_hui$Sigma,
                         TVCL = mp_hui$TVCL, TVV1 = mp_hui$TVV1,
                         TVQ  = mp_hui$TVQ,  TVV2 = mp_hui$TVV2,
                         theta_bsa_cl  = mp_hui$theta_bsa_cl,
                         theta_clcr_cl = mp_hui$theta_clcr_cl,
                         theta_bsa_v1  = mp_hui$theta_bsa_v1,
                         theta_age_q   = mp_hui$theta_age_q,
                         med_bsa    = mp_hui$med_bsa,
                         med_clcr   = mp_hui$med_clcr,
                         med_age    = mp_hui$med_age,
                         height_cm  = mp_hui$height_cm,
                         age_months = mp_hui$age_months) {
  
  sigprop <- as.numeric(Sigma[1, 1])
  sigadd  <- as.numeric(Sigma[2, 2])
  eta     <- as.list(eta)
  eta_m   <- unlist(eta) %>% matrix(nrow = 1)
  cat("hui eta:", paste(names(eta), round(unlist(eta), 4), sep="=", collapse=", "), "\n")
  looppat <- simulate_hui(
    dat_subj,
    TVCL = TVCL, TVV1 = TVV1, TVQ = TVQ, TVV2 = TVV2,
    theta_bsa_cl  = theta_bsa_cl,
    theta_clcr_cl = theta_clcr_cl,
    theta_bsa_v1  = theta_bsa_v1,
    theta_age_q   = theta_age_q,
    med_bsa    = med_bsa,
    med_clcr   = med_clcr,
    med_age    = med_age,
    height_cm  = height_cm,
    age_months = age_months,
    ETAcl  = eta[[1]],
    ETAv2  = eta[[2]]
  )
  
  newobs <- looppat$Conc[which(looppat$time %in% yt)]
  
  sig2j  <- (newobs * sqrt(sigprop) + sqrt(sigadd))^2
  sqwres <- log(sig2j) + (1 / sig2j) * (y - newobs)^2
  
  nOn <- diag(eta_m %*% solve(Omega) %*% t(eta_m))
  return(sum(sqwres) + nOn)
}

pkprof_est_hui <- function(time_and_covs, mp_hui) {
  
  y  <- time_and_covs$dv[which(!is.na(time_and_covs$dv))[1]]
  yt <- round(time_and_covs$time[which(!is.na(time_and_covs$dv))[1]], 1)
  
  init <- c(ETAcl = 0.05, ETAv2 = 0.05)
  
  dat_subj <- subset(time_and_covs, !is.na(dv) | !is.na(amt))
  
  fit <- newuoa(
    par = init,
    fn  = mapbayes_hui,
    dat_subj      = dat_subj,
    y             = y,
    yt            = yt,
    Omega         = mp_hui$Omega,
    Sigma         = mp_hui$Sigma,
    TVCL          = mp_hui$TVCL,
    TVV1          = mp_hui$TVV1,
    TVQ           = mp_hui$TVQ,
    TVV2          = mp_hui$TVV2,
    theta_bsa_cl  = mp_hui$theta_bsa_cl,
    theta_clcr_cl = mp_hui$theta_clcr_cl,
    theta_bsa_v1  = mp_hui$theta_bsa_v1,
    theta_age_q   = mp_hui$theta_age_q,
    med_bsa       = mp_hui$med_bsa,
    med_clcr      = mp_hui$med_clcr,
    med_age       = mp_hui$med_age,
    height_cm     = mp_hui$height_cm,
    age_months    = mp_hui$age_months
  )
  
  return(fit$par)
}


##########################
######### overall ##########
##########################
combine_schedule <- function(dat, patdat, pats) { #estimated concs #what is pats
  dat$id <- 'Population'
  patdat$id <- 'Individual'
  schedule1 <- rbind(pats, patdat, dat)
  isIP1 <- schedule1$id %in% c("Individual","Population")
  schedule1$size <- as.numeric(isIP1)
  schedule1$Profile <- factor(ifelse(isIP1, schedule1$id, 'Simulated'))
  schedule1
}

#calculates concs
getpatdat <- function(dat_subj, mp, pk_est, mp_hui, pk_est_hui, has_conc = TRUE) {
  out_pop <- simulate_mymodel(dat_subj, mp$THETA1, mp$THETA2, mp$THETA3, mp$THETA4, mp$THETA5, mp$THETA6, 
                              ETA1 = 0, ETA2 = 0, ETA3=0, ETA4=0, ETA5 = 0, ETA6 = 0, 
                              beta_scr_cl = mp$beta_scr_cl, beta_bsa_cl = mp$beta_bsa_cl, beta_gender= mp$beta_gender) %>% 
    dplyr::select(-c(A1,A2, A3)) 
  
  out_pop_lower <- simulate_mymodel(dat_subj, THETA1 = 8.72097, THETA2 = 28.0139, THETA3 = 0.732417, THETA4 = 5.08519, THETA5 = 0.109888, THETA6 = 6.88222, 
                                    ETA1 = 0, ETA2 = 0, ETA3=0, ETA4=0, ETA5 = 0, ETA6 = 0, 
                                    beta_scr_cl = 0.416081, beta_bsa_cl = 0.392739, beta_gender= -0.083572) %>% 
    dplyr::select(-c(A1,A2, A3)) 
  
  out_pop_upper <- simulate_mymodel(dat_subj, THETA1 = 10.1582, THETA2 = 33.5705, THETA3 = 1.0256, THETA4 = 6.20806, THETA5 = 0.166228, THETA6 = 16.3414, 
                                    ETA1 = 0, ETA2 = 0, ETA3=0, ETA4=0, ETA5 = 0, ETA6 = 0, 
                                    beta_scr_cl = 0.74403, beta_bsa_cl = 0.945592, beta_gender= 0.110057) %>% 
    dplyr::select(-c(A1,A2, A3)) 
  
  # --- Hui et al. population prediction ---
  out_pop_hui <- simulate_hui(
    dat_subj,
    TVCL = mp_hui$TVCL, TVV1 = mp_hui$TVV1, TVQ = mp_hui$TVQ, TVV2 = mp_hui$TVV2,
    theta_bsa_cl  = mp_hui$theta_bsa_cl,
    theta_clcr_cl = mp_hui$theta_clcr_cl,
    theta_bsa_v1  = mp_hui$theta_bsa_v1,
    theta_age_q   = mp_hui$theta_age_q,
    med_bsa    = mp_hui$med_bsa,
    med_clcr   = mp_hui$med_clcr,
    med_age    = mp_hui$med_age,
    height_cm  = mp_hui$height_cm,
    age_months = mp_hui$age_months,
    ETAcl = 0, ETAv2 = 0
  ) %>% dplyr::select(-c(A1, A2))
  
  out_pop$tag       <- "Population level: Blackman et al."
  dat_subj$tag      <- "Observed Drug Level"
  out_pop_lower$tag <- "Upper 95% population probability: Blackman et al."
  out_pop_upper$tag <- "Lower 95% population probability: Blackman et al."
  out_pop_hui$tag   <- "Population level: Hui et al."
  
  # Harmonise columns before rbind - Blackman (3-cpt) and Hui (2-cpt) have different column sets
  shared_cols <- c("time", "subject_id", "SCR_mmol", "bsa", "pt_gender", "Conc", "evid", "tag")
  align_cols <- function(df) {
    for (col in shared_cols) { if (!col %in% names(df)) df[[col]] <- NA }
    df <- df[, shared_cols]
    # Flatten any list-columns produced by rowwise() %>% mutate() to plain vectors
    for (col in shared_cols) {
      if (is.list(df[[col]])) df[[col]] <- as.numeric(unlist(df[[col]]))
    }
    df
  }
  
  dfs_to_bind <- list(
    align_cols(out_pop),
    align_cols(out_pop_lower),
    align_cols(out_pop_upper),
    align_cols(out_pop_hui)
  )
  
  if (has_conc) {
    IPRED <- simulate_mymodel(dat_subj, mp$THETA1, mp$THETA2, mp$THETA3, mp$THETA4, mp$THETA5, mp$THETA6, 
                              ETA1 = pk_est[1], ETA2 = pk_est[1], 
                              ETA3=pk_est[1], ETA4=pk_est[1], 
                              ETA5 =pk_est[1], ETA6 = pk_est[1], 
                              beta_scr_cl = mp$beta_scr_cl, beta_bsa_cl = mp$beta_bsa_cl, beta_gender= mp$beta_gender) %>% 
      dplyr::select(-c(A1,A2, A3))
    
    IPRED_hui <- simulate_hui(
      dat_subj,
      TVCL = mp_hui$TVCL, TVV1 = mp_hui$TVV1, TVQ = mp_hui$TVQ, TVV2 = mp_hui$TVV2,
      theta_bsa_cl  = mp_hui$theta_bsa_cl,
      theta_clcr_cl = mp_hui$theta_clcr_cl,
      theta_bsa_v1  = mp_hui$theta_bsa_v1,
      theta_age_q   = mp_hui$theta_age_q,
      med_bsa    = mp_hui$med_bsa,
      med_clcr   = mp_hui$med_clcr,
      med_age    = mp_hui$med_age,
      height_cm  = mp_hui$height_cm,
      age_months = mp_hui$age_months,
      ETAcl  = pk_est_hui[1],
      ETAv2  = pk_est_hui[2]
    ) %>% dplyr::select(-c(A1, A2))
    
    IPRED$tag     <- "Individual level: Blackman et al."
    IPRED_hui$tag <- "Individual level: Hui et al."
    
    dfs_to_bind <- c(dfs_to_bind, list(align_cols(IPRED), align_cols(IPRED_hui)))
  }
  
  plt <- do.call(rbind, dfs_to_bind)
  
  p <- ggplot(data=plt, aes(x=time,y=Conc, color=tag)) + geom_line() + 
    geom_point(data=subset(dat_subj, is.na(amt)), aes(x=time, y=dv)) 
  return(list(plt, p))
}

##########################
######### set up model ##########
##########################
setupModel <- function(in_infusmat1 = in_infusmat1, in_concmat = in_concmat, drug, bsa, scr_vals, scr_times, pt_gender, height_cm, dob, has_conc = FALSE) {
  # Compute age in months from date of birth
  age_months <- as.numeric(difftime(Sys.Date(), as.Date(dob), units = "days")) / 30.4375
  
  infN <- nrow(in_infusmat1) #number of infs 1
  conN <- nrow(in_concmat)
  allN <- max(infN, conN)
  
  # check for consistent size and columns
  # check for required columns
  df <- in_infusmat1 #infusion in list. df will have dose time and dur
  ndf <- names(df)
  if(!('dose' %in% ndf)) {
    stop('infusion data should contain "dose" column')
  }
  if(!('time' %in% ndf)) {
    dtvar <- which(vapply(df, inherits, logical(1), 'POSIXt'))[1]
    if(is.na(dtvar)) stop('infusion data should contain date-time variable or "time" column')
    df[,'time'] <- unclass(df[,dtvar]) / 3600 # convert from seconds to hours
  }
  if(!('duration' %in% ndf)) { #sets a default for duration
    df[,'duration'] <- 2 # default: 2 hours
  }
  df[is.na(df[,'duration']),'duration'] <- 2
  df <- df[,c('dose','time','duration')]
  in_infusmat1 <- df[complete.cases(df),]
  
  df <- in_concmat
  ndf <- names(df)
  if(!('dose' %in% ndf)) {
    stop('concentration data should contain "dose" column')
  }
  if(!('time' %in% ndf)) {
    dtvar <- which(vapply(df, inherits, logical(1), 'POSIXt'))[1]
    if(is.na(dtvar)) stop('concentration data should contain date-time variable or "time" column')
    df[,'time'] <- unclass(df[,dtvar]) / 3600 # convert from seconds to hours
  }
  df <- df[,c('dose','time')]
  in_concmat <- df[complete.cases(df),]
  
  # Convert all times to hours first
  in_infusmat1$time <- as.numeric(in_infusmat1$time) / 3600
  in_concmat$time   <- as.numeric(in_concmat$time)   / 3600
  scr_times         <- as.numeric(scr_times)         / 3600
  
  # has_conc is passed in from the app (detected before checkInputMatrix adds fallback default)
  # Then compute time0 and zero everything.
  # Only include in_concmat times when real conc data exists; the default fallback
  # timestamp is near epoch-now and causes time0 to be wildly wrong otherwise.
  time0 <- if (has_conc) {
    min(in_infusmat1$time, in_concmat$time)
  } else {
    min(in_infusmat1$time)
  }
  in_infusmat1$time <- in_infusmat1$time - time0
  in_concmat$time   <- in_concmat$time   - time0
  scr_times = round(as.numeric(scr_times - time0),1)
  # data check for SCR values
  if(is.na(scr_times[1])){
    scr_times[1] = 0
  }
  # Drop entries where either the time or the value is missing
  valid_scr <- !is.na(scr_times) & !is.na(scr_vals)
  scr_vals  = scr_vals[valid_scr]
  scr_times = scr_times[valid_scr]
  last.t <- if (has_conc) {
    max(in_infusmat1$time, in_concmat$time, scr_times, na.rm=TRUE)
  } else {
    max(in_infusmat1$time, scr_times, na.rm=TRUE)
  }
  usrendt <- max(last.t + 12, 72) #max for input
  # check covariates
  negNA <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x[x < 0] <- NA
    x
  }
  usrgender <- usrSCR <- usrbsa <- NA
  
  usrgender  <- pt_gender
  usrSCR     <- scr_vals * 88.4   # mg/dL -> mmol/L
  usrbsa     <- bsa
  mp         <- model_mtx(usrbsa = usrbsa, usrSCR = usrSCR, usrgender = usrgender)
  mp_hui     <- model_hui(usrbsa = usrbsa, usrSCR_mgdl = scr_vals, usr_height_cm = height_cm, usr_age_months = age_months)
  
  # construct dataframe
  timevec <- seq(from = 0, to = usrendt, by = 0.1)
  print(usrendt)
  time_and_covs <- data.frame(time = timevec, amt = NA, dv = NA, rate = NA,
                              mdv = 1, evid = 0, dur = NA,
                              pt_gender = mp$pt_gender, bsa = mp$bsa,
                              SCR_mmol = NA, SCR_mgdl = NA)
  
  # Replace the SCR fill loop (lines 653-658):
  for (i in 1:length(scr_times)){
    if (is.na(scr_times[i]) || is.na(usrSCR[i])) next
    # snap to nearest timevec row instead of exact match
    diffs <- abs(time_and_covs$time - scr_times[i])
    row <- which(diffs == min(diffs))[1]
    if (length(row) > 0) {
      time_and_covs$SCR_mmol[row] <- usrSCR[i]
      time_and_covs$SCR_mgdl[row] <- usrSCR[i] / 88.4
    }
  }
  
  
  inf = in_infusmat1
  for (i in 1:nrow(inf)){
    time = inf$time[i]
    dose = inf$dose[i]
    duration = inf$duration[i]
    rate = dose/duration
    idx = which(abs(time_and_covs$time - round(time,1)) < 0.01)
    time_and_covs$amt[idx] = dose ; time_and_covs$dur[idx] = duration ; time_and_covs$rate[idx] = rate
    time_and_covs$evid[idx] = 1
  }
  conc = in_concmat
  if (has_conc) {
    for (i in 1:nrow(conc)){
      row = i
      time = conc$time[i]
      dv = conc$dose[i]
      idx = which(abs(time_and_covs$time - round(time,1)) < 0.01)
      time_and_covs$dv[idx] = dv ; time_and_covs$mdv[idx] = 0
    }
  }
  
  library(tidyverse)
  time_and_covs <- time_and_covs %>% 
    fill(dur, SCR_mmol, SCR_mgdl, .direction = "downup")
  
  dat_subj = subset(time_and_covs, !is.na(dv) | !is.na(amt))
  
  
  def.par <- list(
    y = conc[,'dose'], yt = conc[,'time']
  )
  
  # pk_est <- pkprof_est_mymodel(time_and_covs, mp)
  # pk_est_hui <- pkprof_est_hui(time_and_covs, mp_hui = mp_hui)
  
  # With this:
  if (has_conc) {
    pk_est     <- pkprof_est_mymodel(time_and_covs, mp)
    print("pk_est my model ran")
    
    pk_est_hui <- pkprof_est_hui(time_and_covs, mp_hui = mp_hui)
    print("pk_est Hui et al. ran")
    
  } else {
    pk_est     <- c(ETAcl=0, ETAv1=0, ETAq2=0, ETAv2=0, ETAq3=0, ETAv3=0)
    print("pk_est my model ran")
    
    pk_est_hui <- c(ETAcl=0, ETAv2=0)
    print("pk_est Hui et al. ran")
    
  }
  
  preds <- getpatdat(dat_subj = time_and_covs, mp, pk_est, mp_hui, pk_est_hui, has_conc = has_conc)
  
  concdat <- if (has_conc) {
    data.frame(y = in_concmat[,'dose'], yt = in_concmat[,'time'], colour = " Observed Drug Level")
  } else {
    data.frame(y = numeric(0), yt = numeric(0), colour = character(0))
  }
  
  return(list(concdat, preds[[1]], preds[[2]], mp))
}

makePlots <- function(schedule1, concdat) {
  level_order <- c(
    "Observed Drug Level",              # Col 1, Row 1
    " ",                                # Col 2, Row 1 (The Blank)
    "Serum Creatinine (mg/dL)",          # Col 1, Row 2
    "Glucarpidase Consensus Guidelines", # Col 2, Row 2
    "Population level: Blackman et al.",
    "Individual level: Blackman et al.",
    "Population level: Hui et al.",
    "Individual level: Hui et al.",
    "Upper 95% population probability: Blackman et al.",
    "Lower 95% population probability: Blackman et al."
  )
  
  cols <- c(
    "Observed Drug Level" = "red",
    " " = "#00000000",
    "Serum Creatinine (mg/dL)" = "purple",
    "Glucarpidase Consensus Guidelines" = "blue2",
    "Population level: Blackman et al." = "deeppink4",
    "Individual level: Blackman et al." = "deeppink1",
    "Population level: Hui et al." = "darkorange1",
    "Individual level: Hui et al." = "darkgoldenrod2",
    "Upper 95% population probability: Blackman et al." = "grey50",
    "Lower 95% population probability: Blackman et al." = "grey50"
  )
  
  # Prep datasets with factors to enforce order in Plotly
  gluc_guidelines <- data.frame(
    time = c(24, 36, 42, 48), 
    amt = c(50, 30, 10, 5) * 0.454,
    tag = factor("Glucarpidase Consensus Guidelines", levels = level_order)
  ) 
  
  df_ribbon <- schedule1 %>%
    filter(tag %in% c("Lower 95% population probability: Blackman et al.", 
                      "Upper 95% population probability: Blackman et al.")) %>%
    dplyr::select(time, tag, Conc) %>% 
    pivot_wider(names_from = tag, values_from = Conc) %>% 
   # rename(lower = 2, upper = 3)
    rename(lower = "Lower 95% population probability: Blackman et al.",
           upper = "Upper 95% population probability: Blackman et al.")
  mult <- 30
  scr_plot_data <- schedule1 %>%
    dplyr::distinct(time, SCR_mmol, .keep_all = TRUE) %>%
    dplyr::arrange(time) %>%
    filter(!is.na(SCR_mmol)) %>%
    filter(SCR_mmol != lag(SCR_mmol, default = -Inf)) %>%
    mutate(SCR_mgdl = SCR_mmol / 88.4,
           y_scaled = SCR_mgdl * mult,
           tag = factor("Serum Creatinine (mg/dL)", levels = level_order))
  
  # Ensure the main data is factored
  schedule1$tag <- factor(schedule1$tag, levels = level_order)
  
  ggplot() +
    # Ribbon first
    geom_ribbon(data = df_ribbon, aes(x = time, ymin = lower, ymax = upper), 
                fill = "grey10", alpha = 0.15) +
    
    # Lines - use 'name = tag' for Plotly
    geom_line(data = schedule1 %>% filter(tag %in% level_order[4:10]), 
              aes(x = time, y = Conc, color = tag, group = tag, name = tag)) +
    
    # Observed Points
    geom_point(data = concdat, 
               aes(x = yt, y = y, color = factor("Observed Drug Level", levels = level_order), 
                   name = factor("Observed Drug Level", levels = level_order)), 
               size = 2.5) +
    
    # Glucarpidase Points
    geom_point(data = gluc_guidelines, aes(x = time, y = amt, color = tag, name = tag), 
               shape = 5, size = 2.5, stroke = 0.6) +
    
    # Serum Creatinine Points
    geom_point(data = scr_plot_data, aes(x = time, y = y_scaled, color = tag, name = tag), 
               size = 3, shape = 2, stroke = 0.6) +
    
    # FIX: Removed alpha = 0. The color vector will handle transparency via "#00000000"
    geom_point(data = data.frame(x = 0, y = 0.1), 
               aes(x = x, y = y, color = factor(" ", levels = level_order)), 
               show.legend = TRUE) +
    
    scale_colour_manual(name = NULL, values = cols, breaks = level_order, drop = FALSE) +
    scale_x_continuous(name = "Time (hours)", breaks = seq(0, 72, by = 12)) +
    scale_y_continuous(name = "Concentration (mg/L)")
}



arrange_data <- function(df) {
  dfn <- tolower(names(df))
  rateCol <- grep('rate', dfn)
  amntCol <- grep('amt', dfn)
  concCol <- grep('conc', dfn)
  duraCol <- grep('dur', dfn)
  timeCol <- grep('time', dfn)
  timeVar <- 'time'
  if(length(rateCol) == 0) stop('rate column should exist')
  if(length(amntCol) == 0) stop('amt column should exist')
  if(length(concCol) == 0) stop('conc column should exist')
  if(length(timeCol) == 0) {
    addlCol <- setdiff(seq(ncol(df)), c(rateCol, amntCol, concCol, duraCol, timeCol))
    dtCol <- character(0)
    # is additional column a date-time variable?
    if(length(addlCol) > 0) {
      isDT <- apply(df[,addlCol,drop=FALSE], 2, function(i) tryCatch(pkdata::guessDateFormat(i), error=function(e) NA_character_))
      dtCol <- dfn[addlCol[!is.na(isDT)]]
    }
    if(length(dtCol) != 1) stop('date-time or time column should exist')
    timeVar <- 'datetime'
    df[,timeVar] <- as.POSIXct(df[,dtCol], format = isDT[!is.na(isDT)])
    timeCol <- match(timeVar, names(df))
  }
  df_i <- df[!is.na(df[,rateCol]),c(rateCol,timeCol,duraCol)]
  df_c <- df[!is.na(df[,concCol]),c(concCol,timeCol)]
  names(df_i)[1:2] <- c('dose',timeVar)
  if(ncol(df_i) == 3) names(df_i)[3] <- 'duration'
  names(df_c) <- c('dose',timeVar)
  covar <- list(
    wt = df[1,grep('weight', dfn)],
    alb = df[1,grep('alb', dfn)],
    ada = df[1,grep('ada', dfn)]
  )
  list(infList = list(df_i), conList = list(df_c), covariates = covar)
}