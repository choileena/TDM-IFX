co <- function(alpha,beta,D,k,t) {
#   D*((alpha-k)/(alpha-beta))*exp(-alpha*t) + D*((k-beta)/(alpha-beta))*exp(-beta*t)
  ((alpha-k)/(alpha-beta))*exp(-alpha*t) + ((k-beta)/(alpha-beta))*exp(-beta*t)
}
# coinf <- function(alpha,beta,D,t,A,B) D*(A/alpha*(1-exp(-alpha*t)) + B/beta*(1-exp(-beta*t)) )
coinf <- function(alpha,beta,D,t,A,B) (A/alpha*(1-exp(-alpha*t)) + B/beta*(1-exp(-beta*t)) )

nwn <- function(x) suppressWarnings(as.numeric(x))

# checkInputMatrix <- function(m, makeNumeric = TRUE) {
#   #check rows for any empty ("") values or characters
#   badrows <- c()
#   for(i in 1:nrow(m)){
#     suppressWarnings(g <- any(is.na(as.numeric(m[i,]))))
#     if(g) badrows <- c(badrows,i)
#   }
#   if(length(badrows>=1)) m <- m[-badrows,]
#   if(is.vector(m)) m <- matrix(m, nrow=1) 
#   if(nrow(m)==0) m <-  matrix(c(0,0,0), nrow=1) 
#   m <- matrix(as.numeric(m), nrow=nrow(m))
#   m <- m[rowSums(nchar(m) == 0) == 0,,drop = FALSE]
#   if(makeNumeric) mode(m) <- 'numeric'
#   
#   #check for any negative numbers
#   m[m<0] <- 0
#   
#   m
# }
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
  m[,'time'] <- nwn(m[,'dt']) / 3600
  m <- m[complete.cases(m),]
  if(nrow(m) == 0) {
    m <- def
  }
  m
}

getPatDat <- function(clvec=1, v1vec=1, qvec=0, v2vec=1,
                      usrendt, usrbdose = NULL, timevec, usrbt = NULL, usrwt,
                      usriet, usrist, usridose, Omega){
  # usrendt: user end date-time
  # usrbdose: bolus dose
  # timevec: time sequence (redundant)
  # usrbt: bolus time
  # usrwt: user weight
  # usriet: infusion end time
  # usrist: infusion start time
  # usridose: infusion dose
  timevec = seq(0, usrendt + (2*24*7), by=1)
  ltv <- length(timevec)
  pats2 <- vector('list', ltv)
  dur <- usriet - usrist
  dosewgtByDur <- usridose * usrwt / dur
  begt_j <- match(usrist, timevec)
  endt_j <- match(usriet, timevec)
  for(i in 1:length(clvec)) {
    concvec2 <- rep(0,ltv)
    k10 <- clvec[i] / v1vec[i]
    k12 <- qvec[i] / v1vec[i]
    k21 <- qvec[i] / v2vec[i]
    if(nrow(Omega)==2) {k12=k21=0}
    a <- .5*(k12 + k21 + k10 + sqrt( (k12+k21+k10)^2 -4*k21*k10 ))
    b <- k12 + k21 + k10-a
    A <- (k21-a)/(b-a)
    B <- (k21-b)/(a-b)
    if(nrow(Omega)==2) {b=b+.001}
    coinf_ex <- coinf(alpha=a, beta=b, D=NA, t=timevec, A=A, B=B)
    tva <- exp(-a * timevec)
    tvb <- exp(-b * timevec)
    calc_decay <- function(w, d) w * (A/a * (1-exp(-a*d)) * tva + B/b * (1-exp(-b*d)) * tvb )
    decay_ij <- mapply(calc_decay, dosewgtByDur, dur)
    for(j in 1:length(usridose)){
      infvec <- rep(0,ltv)
      infusion <- coinf_ex * dosewgtByDur[j]
      begt <- begt_j[j]
      endt <- endt_j[j]
      infvec[begt:endt] <- infusion[seq(endt-begt+1)]

      if(endt != ltv){
        decay <- decay_ij[,j]
        infvec[(endt+1):ltv] <- decay[2:(ltv-endt+1)]
      }
      concvec2 <- concvec2 + infvec / v1vec[i]
    }
    pats2[[i]] <- data.frame(id=i, t=timevec, concentration = concvec2)
  }
  do.call(rbind, pats2)
}

mapbayes <- function(eta,
                     y,
                     yt,
                     Omega,
                     Sigma,
                     model,
                     ncpt,
                     usrendt,timevec,usrwt,
                     usriet,usrist,usridose) {
  sigprop <- as.numeric(Sigma[1,1])
  sigadd <- as.numeric(Sigma[2,2])
  clloop <- model[1] * exp(eta['ETAcl'])
  if(ncpt == 2) {
    v1loop <- model[2] * exp(eta['ETAv1'])
    qloop <- model[3] * exp(eta['ETAq'])
    v2loop <- model[4] * exp(eta['ETAv2'])
  } else {
    v1loop <- model[2] * exp(eta['ETAv'])
    qloop <- 0
    v2loop <- 1
  }

  looppat <- getPatDat(clvec=clloop,v1vec=v1loop,qvec=qloop,v2vec=v2loop,
                       usrendt=usrendt,usrbdose=NULL,timevec=NA,usrbt=NULL,usrwt=usrwt,
                       usriet=usriet,usrist=usrist,usridose=usridose,Omega=Omega)
  newobs = looppat$concentration[looppat$t %in% yt]

  # http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3339294/
  sig2j <- (newobs*sqrt(sigprop) + sqrt(sigadd))^2    #<------------------Check against pubmed
  sqwres <- log(sig2j) + (1/sig2j)*(y-newobs)^2

  eta_m <- matrix(eta, nrow=1)
  nOn <- diag(eta_m %*% solve(Omega) %*% t(eta_m))
  return(sum(sqwres) + nOn)
}

bldplot <- function(dat, title, tdrupper, tdrlower) {
  mx <- max(dat[,'t'], na.rm = TRUE)  
  xby <- 1 # (24 * 7) # HERE (24 * 7)
  # set y-limit for concentration; exclude 0 as log10(0) = -Inf
  min_conc <- min(dat[dat$concentration > 1e-6, 'concentration'])
  ggplot(data = dat) + ggtitle(title) +
    geom_line(aes(x=t/(24 * 7), y=concentration, color=Profile, group=as.factor(id), size=as.factor(size))) + # HERE x=t/(24 * 7), check size
    scale_size_manual(values=c(.1,1)) +
    scale_y_log10(limits = c(min_conc, NA)) +   # HERE log scale, set limit
    scale_x_continuous(breaks = seq(xby, mx, by = xby)) +
    geom_hline(yintercept=tdrupper, linetype="dashed") +
    geom_hline(yintercept=tdrlower, linetype="dashed") +
    labs(y = 'Infliximab conc. (mcg/mL)', # HERE
         x = 'Time (weeks)',
         title = '') + 
    # coord_cartesian(ylim=c(0.01, 100)) + 
    # expand_limits(y=0.1) +
    guides(group = "none", size = "none") + xlim(0, 16)
}

pkprof_est <- function(model_params) {
  ncpt <- model_params$ncpt
  # removed 'usrbdose' & 'usrbt'
  mp <- model_params[c('y','yt','Omega','Sigma','ncpt','usrendt','timevec','usrwt','usriet','usrist','usridose')]
  mp$fn <- mapbayes
  if(ncpt == 2) {
    init <- c(ETAcl=.05, ETAv1=.05, ETAq=.05, ETAv2=.05)
    cl0 <- model_params$cl0
    v10 <- model_params$v10
    q0 <- model_params$q0
    v20 <- model_params$v20
    mp$model <- c(cl0, v10, q0, v20)
  } else if(ncpt == 1) {
    init <- c(ETAcl=.05, ETAv=.05)
    cl0 <- model_params$cl0
    v10 <- model_params$v0
    q0 <- 0
    v20 <- 1
    mp$model <- c(cl0, v10)
  }
  mp$par <- init
  fit <- do.call(minqa::newuoa, mp)
  exp_par <- exp(fit$par)
  if(ncpt == 1) {
    exp_par[3:4] <- c(0,1)
  }
  c(c(cl0, v10, q0, v20) * exp_par[1:4], c(cl0, v10, q0, v20))
}

model_infliximab <- function(usrN, usrwt, usralb, usrada) {
  clpop = 0.296/(24) # = 0.0123 ; 0.296 L/day 
  v1pop = 3.3
  qpop = 0.0781/(24) # =0.003254 ; 0.0781 L/day 
  v2pop = 1.16
  OmegaCL = 0.313
  OmegaV1 = 0.0985
  OmegaQ = 1.11
  OmegaV2 = 0.761
  properr = 0.419
  adderr = 0
  b_wt_cl = 0.612
  b_wt_v1 = 0.696
  b_wt_q = 1.15
  b_wt_v2 = 0.604
  beta_alb_cl = -2.3 
  beta_ada_cl = 0.231
  Omega = diag(c(OmegaCL,OmegaV1,OmegaQ,OmegaV2))
  Sigma = diag(c(properr,adderr))
  ETA_mat = phonTools::rmvtnorm(n=usrN,k=nrow(Omega),means=rep(0,nrow(Omega)),sigma=Omega^2)
  scale_wgt <- usrwt / 70
  scale_alb <- usralb / 4
  cl0 = clpop * scale_wgt^b_wt_cl * (scale_alb)^beta_alb_cl * (1 + usrada * beta_ada_cl) 
  v10 = v1pop * scale_wgt^b_wt_v1 
  q0 = qpop * scale_wgt^b_wt_q 
  v20 = v2pop * scale_wgt^b_wt_v2
  clvec = cl0 * exp(ETA_mat[,1])
  v1vec = v1pop * scale_wgt^b_wt_v1 * exp(ETA_mat[,2])
  qvec = qpop * scale_wgt^b_wt_q * exp(ETA_mat[,3])
  v2vec = v2pop * scale_wgt^b_wt_v2 * exp(ETA_mat[,4])
  list(ncpt = 2, Omega = Omega, Sigma = Sigma, tdrupper = 10, tdrlower = 5, usrwt = usrwt,
    clvec = clvec, v1vec = v1vec, qvec = qvec, v2vec = v2vec, cl0 = cl0, v10 = v10, q0 = q0, v20 = v20
  )
}

combine_schedule <- function(dat, patdat, pats) {
    dat$id <- 'Population'
    patdat$id <- 'Individual'
    schedule1 <- rbind(pats, patdat, dat)
    isIP1 <- schedule1$id %in% c("Individual","Population")
    schedule1$size <- as.numeric(isIP1)
    schedule1$Profile <- factor(ifelse(isIP1, schedule1$id, 'Simulated'))
    schedule1
}

mk_patdat_bidtid <- function(pkprof) {
    #simulate 100 days for induction phase
    nDays <- 100 # 180
    day10 <- nDays + 10
    # inf.st.time4 <- seq(0, nDays*24, by = 7 * 4 * 24)
    # inf.st.time6 <- seq(0, nDays*24, by = 7 * 6 * 24)
    # inf.st.time8 <- seq(0, nDays*24, by = 7 * 8 * 24)

    inf.st.time.ind <- c(0, 2*7*24, 6*7*24, 14*7*24) # HERE LC 
    # len4 <- length(inf.st.time)
    # len6 <- length(inf.st.time6)
    # len8 <- length(inf.st.time8)

    len.ind <- length(inf.st.time.ind) # HERE LC 
    
    def.args <- list(clvec=pkprof$p[1], v1vec=pkprof$p[2], qvec=pkprof$p[3], v2vec=pkprof$p[4],
      usrendt=day10*24, usrbdose=NULL, timevec=NULL,
      usrbt=NULL, usrwt=pkprof$wt, Omega=pkprof$Omega
    )
    # "usriet" is end time; assumed duration of 2 hours
    # arg4 <- list(usriet=inf.st.time4+2, usrist=inf.st.time4)
    # arg6 <- list(usriet=inf.st.time6+2, usrist=inf.st.time6)
    # arg8 <- list(usriet=inf.st.time8+2, usrist=inf.st.time8)

    arg.ind <- list(usriet=inf.st.time.ind+2, usrist=inf.st.time.ind) # HERE LC 

    # r5q4 <- do.call(getPatDat, c(def.args, arg4, list(usridose=rep(5,len4))))
    # r5q6 <- do.call(getPatDat, c(def.args, arg6, list(usridose=rep(5,len6))))
    # r5q8 <- do.call(getPatDat, c(def.args, arg8, list(usridose=rep(5,len8))))
    # r75q4 <- do.call(getPatDat, c(def.args, arg4, list(usridose=rep(7.5,len4))))
    # r75q6 <- do.call(getPatDat, c(def.args, arg6, list(usridose=rep(7.5,len6))))
    # r75q8 <- do.call(getPatDat, c(def.args, arg8, list(usridose=rep(7.5,len8))))
    # r10q4 <- do.call(getPatDat, c(def.args, arg4, list(usridose=rep(10,len4))))
    # r10q6 <- do.call(getPatDat, c(def.args, arg6, list(usridose=rep(10,len6))))
    # r10q8 <- do.call(getPatDat, c(def.args, arg8, list(usridose=rep(10,len8))))

    r5.ind <- do.call(getPatDat, c(def.args, arg.ind, list(usridose=rep(5,len.ind)))) # HERE LC 
    r75.ind <- do.call(getPatDat, c(def.args, arg.ind, list(usridose=rep(7.5,len.ind)))) # HERE LC 
    r10.ind <- do.call(getPatDat, c(def.args, arg.ind, list(usridose=rep(10,len.ind)))) # HERE LC 

    # list(r5q4, r75q4, r10q4, r5q6, r75q6, r10q6, r5q8, r75q8, r10q8)
    
    list(r5.ind, r75.ind, r10.ind) # HERE LC 
}

getConcRange <- function(dat, bt, tt) {
  tmpt <- dat[,'t']
  ix <- which(tmpt >= bt & tmpt < tt)
  rev(range(dat[ix,'concentration']))
}

setupModel <- function(infList = NULL, conList = NULL, drug, ...) {
  if(inherits(infList, 'data.frame')) {
    infList <- list(infList)
  }
  if(inherits(conList, 'data.frame')) {
    conList <- list(conList)
  }
  infN <- length(infList)
  conN <- length(conList)
  allN <- max(infN, conN)

  # check for consistent size and columns
  if(infN > 1) {
    chk <- vapply(infList, function(i) length(setdiff(names(i), names(infList[[1]]))), numeric(1))
    if(sum(chk) > 0) stop('all data.frames in infList should be compatible')
    if(infN < allN) stop('not enough data.frames provided in infList')
  }
  if(conN > 1) {
    chk <- vapply(conList, function(i) length(setdiff(names(i), names(conList[[1]]))), numeric(1))
    if(sum(chk) > 0) stop('all data.frames in conList should be compatible')
    if(conN < allN) stop('not enough data.frames provided in conList')
  }

  # check for required columns
  for(i in seq_len(infN)) {
    df <- infList[[i]]
    ndf <- names(df)
    if(!('dose' %in% ndf)) {
      stop('infusion data should contain "dose" column')
    }
    if(!('time' %in% ndf)) {
      dtvar <- which(vapply(df, inherits, logical(1), 'POSIXt'))[1]
      if(is.na(dtvar)) stop('infusion data should contain date-time variable or "time" column')
      df[,'time'] <- unclass(df[,dtvar]) / 3600 # convert from seconds to hours
    }
    if(!('duration' %in% ndf)) {
      df[,'duration'] <- 2 # default: 2 hours
    }
    df[is.na(df[,'duration']),'duration'] <- 2
    df <- df[,c('dose','time','duration')]
    infList[[i]] <- df[complete.cases(df),]
  }
  for(i in seq_len(conN)) {
    df <- conList[[i]]
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
    conList[[i]] <- df[complete.cases(df),]
  }

  # determine MIN time and fix offset
  infMin <- suppressWarnings(min(vapply(infList, function(i) min(i[,'time']), numeric(1))))
  conMin <- suppressWarnings(min(vapply(conList, function(i) min(i[,'time']), numeric(1))))
  time0 <- min(infMin, conMin)
  for(i in seq_len(infN)) {
    inftime <- infList[[i]][,'time'] - time0
    infList[[i]][,'time'] <- inftime
    infList[[i]][,'endtime'] <- inftime + infList[[i]][,'duration']
  }
  for(i in seq_len(conN)) {
    conList[[i]][,'time'] <- conList[[i]][,'time'] - time0
  }

  # determine MAX time
  infMax <- suppressWarnings(max(vapply(infList, function(i) max(i[,'endtime']), numeric(1))))
  conMax <- suppressWarnings(max(vapply(conList, function(i) max(i[,'time']), numeric(1))))
  last.t <- max(infMax, conMax)
  usrendt <- last.t + 12

  # check covariates
  negNA <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x[x < 0] <- NA
    x
  }
  covariates <- list(...)
  ncov <- names(covariates)
  usrwt <- usralb <- usrada <- NA
  if('wt' %in% ncov) {
    usrwt <- negNA(covariates[['wt']])
  }
  if('alb' %in% ncov) {
    usralb <- negNA(covariates[['alb']])
  }
  if('ada' %in% ncov) {
    usrada <- negNA(covariates[['ada']])
  }

  usrN <- 10
  if(drug=="infliximab"){
    mp <- model_infliximab(usrN, usrwt, usralb, usrada)
  }
  Omega <- mp$Omega
  timevec <- NULL
  def.par <- list(
    y = conList[[1]][,'dose'], yt = conList[[1]][,'time'],
    usrendt = usrendt, timevec = timevec,
    usriet = infList[[1]][,'endtime'], usrist = infList[[1]][,'time'], usridose = infList[[1]][,'dose']
  )
  pk_est <- pkprof_est(c(def.par, mp))
  data_args <- list(clvec=pk_est[5], v1vec=pk_est[6], qvec=pk_est[7], v2vec=pk_est[8])
  pdat_args <- list(clvec=pk_est[1], v1vec=pk_est[2], qvec=pk_est[3], v2vec=pk_est[4])
  pats_args <- list(clvec=mp$clvec, v1vec=mp$v1vec, qvec=mp$qvec, v2vec=mp$v2vec)
  schedule <- vector('list', allN)
  for(i in seq_along(schedule)) {
    s_args <- list(
      usrendt=usrendt, timevec=timevec, usrwt=usrwt, Omega=Omega,
      usriet=infList[[i]][,'endtime'], usrist=infList[[i]][,'time'], usridose=infList[[i]][,'dose']
    )
    dat_schedule1 <- do.call(getPatDat, c(data_args, s_args))
    patdat_schedule1 <- do.call(getPatDat, c(pdat_args, s_args))
    pats_schedule1 <- do.call(getPatDat, c(pats_args, s_args))
    schedule[[i]] <- combine_schedule(dat_schedule1, patdat_schedule1, pats_schedule1)
  }
  tdrupper <- mp$tdrupper
  tdrlower <- mp$tdrlower
  concdat <- data.frame(y = conList[[1]][,'dose'], yt = conList[[1]][,'time'], colour = " Observed Drug Level")
  PKprof <- list(Omega = Omega, p = unname(unlist(pdat_args)), wt = usrwt)
  list(schedule, tdrupper, tdrlower, concdat, PKprof)
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
