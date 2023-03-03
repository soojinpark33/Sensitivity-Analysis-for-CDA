# Population data generation
popDataGen <- function(N = 1000000, effectSize = c("small", "medium", "large"), seed = NULL){
  
  # Set parameters
  p0 <- 0.5       # target prop. of White men (R = 0)
  lb.C <- 25      # lower bound for C
  ub.C <- 75      # upper bound for C
  avg.C <- 50     # truncated normal mean for C
  sd.C <- 12      # truncated normal sd for C
  shift.C <- 2    # difference in truncated normal means between two groups (R = 0, 1) for C
  mean.X1 <- 0.5  # poisson mean of S for group R = 0
  mean.X2 <- 0.8  # poisson mean of S for group R = 1
  err.fitM <- err.fitX <- err.fitU <- 0.5   # for error term of M model
  err.fitY <- 0.5                           # for error term of Y model
  
  # Create variables
  set.seed(seed)
  R <- as.numeric(rbinom(N, size = 1, prob = p0))
  C <- rep(NA, N)
  C[R == 0] <- rtruncnorm(sum(R == 0), a = lb.C, b = ub.C, mean = avg.C, sd = sd.C)
  C[R == 1] <- rtruncnorm(sum(R == 1), a = lb.C, b = ub.C, mean = avg.C - shift.C, sd = sd.C)
  err.M <- rnorm(N, sd = err.fitM)
  err.X <- rnorm(N, sd = err.fitX)
  err.U <- rnorm(N, sd = err.fitU)
  err.Y <- rnorm(N, sd = err.fitY)
  C.grp <- cut(C, breaks = c(0, 50, Inf), include.lowest = TRUE, labels = c(1, 2))
  X <- coef.X[1] + coef.X[2]*R + coef.X[3]*ifelse(C.grp == 2, 1, 0) + err.X 
  U <- coef.U[1] + coef.U[2]*R + coef.U[3]*X + coef.U[4]*ifelse(C.grp == 2, 1, 0) + err.U 
  C.grp.ref2 <- relevel(C.grp, ref = 2)
  M <- coef.M[1] + coef.M[2]*R + coef.M[3]*X + coef.M[4]*U + coef.M[5]*ifelse(C.grp == 2, 1, 0) + err.M
  Y <- coef.Y[1] + coef.Y[2]*R + coef.Y[3]*X + coef.Y[4]*M + coef.Y[5]*U + coef.Y[6]*ifelse(C.grp == 2, 1, 0) + err.Y
  pop.dat <- data.frame(Y, R, X, M, U, C, C.grp, C.grp.ref2)
  names(pop.dat) <- c("Y", "R", "X", "M", "U", "C", "C.grp", "C.grp.ref2")
  pop.dat$R <- as.factor(pop.dat$R)
  
  coef.M0 <- coef.M
  coef.Y0 <- coef.Y
  
  # Adjust data for pre-specified effect 
  if(!is.null(effectSize)){
    if(effectSize == "small"){
      coef.M0[4] <- -0.1418; coef.Y0[5] <- -0.1439
    } else if (effectSize == "medium"){
      coef.M0[4] <- -0.3859; coef.Y0[5] <- -0.4141
    } else if (effectSize == "large"){
      coef.M0[4] <- -0.5926; coef.Y0[5] <- -0.6891
    }
  }
  pop.dat$M <- coef.M[1] + coef.M[2] * R + coef.M[3] * X + coef.M0[4] * U +
    coef.M[5] * ifelse(C.grp == 2, 1, 0) + err.M
  pop.dat$Y <- coef.Y[1] + coef.Y[2] * R + coef.Y[3] * X + coef.Y[4] * M +
    coef.Y0[5] * U + coef.Y[6] * ifelse(C.grp == 2, 1, 0) + err.Y
  rsq.YU <- (anova(lm(pop.dat$Y ~ R + X + pop.dat$M + C.grp))[5, 2] -
               anova(lm(pop.dat$Y ~ U + R + X + pop.dat$M + C.grp))[6, 2])/
    anova(lm(pop.dat$Y ~ R + X + pop.dat$M + C.grp))[5, 2]
  rsq.MU <- (anova(lm(pop.dat$M ~ R + X + C.grp))[4, 2] -
               anova(lm(pop.dat$M ~ U + R + X + C.grp))[5, 2])/
    anova(lm(pop.dat$M ~ R + X + C.grp))[4, 2]
  
  y0 <- mean(pop.dat$Y[pop.dat$R == 0 & pop.dat$C.grp == 2])
  y1 <- mean(pop.dat$Y[pop.dat$R == 1 & pop.dat$C.grp == 2])
  m0 <- pop.dat$M[pop.dat$R == 0 & pop.dat$C.grp == 2]
  m1 <- pop.dat$M[pop.dat$R == 1 & pop.dat$C.grp == 2]
  m01 <- sample(m0, (sum(pop.dat$R == 1 & pop.dat$C.grp == 2)), replace = T)
  ym0 <- mean(coef.Y[1] + coef.Y[2] * 1 + coef.Y[3] * pop.dat$X[pop.dat$R == 1 & pop.dat$C.grp == 2] +
                coef.Y[4] * m01 + coef.Y0[5] * pop.dat$U[pop.dat$R == 1 & pop.dat$C.grp == 2] +
                coef.Y[6] * 1 + err.Y[pop.dat$R == 1 & pop.dat$C.grp == 2])
  Red.true <- y1 - ym0
  Rem.true <- ym0 - y0
  #(y1 - y0)
  
  out <- list(data = pop.dat, ry = rsq.YU, rm = rsq.MU, Red.true = Red.true, Rem.true = Rem.true)
  return(out)
  
}

# CI computation
sens.for.se <- function(boot.res, fit.y, fit.m = NULL, mediators = NULL, covariates, treat, sel.lev.treat,
                        conf.level = 0.95, ry, rm){
  
  # Match arguments
  call <- match.call()
  data <- model.frame(fit.y)
  
  # Warning for inappropriate settings
  if(!(sel.lev.treat %in% levels(data[, treat]))){
    stop("'sel.lev.treat' must be one of the levels of treatment")
  }
  if(sel.lev.treat == levels(data[, treat])[1]){
    stop("'sel.lev.treat' must not be the reference level of treatment")
  }
  
  # Extract outcome model and data from bootstrapping results
  outcome <- all.vars(formula(fit.y))[[1]]
  len.treat <- length(levels(data[, treat]))
  # new data with releveled treatment
  data.new <- data
  data.new[, treat] <- relevel(data.new[, treat], ref = sel.lev.treat)
  # new outcome model with releveled treatment
  fit.y.new <- update(fit.y, data = data.new)
  lev <- which(levels(data[, treat]) == sel.lev.treat)
  
  if(is.null(fit.m) & !is.null(mediators)){
    
    # Possibility of multiple mediators
    if(length(class(fit.m)) == 1 && class(fit.m) == "list"){
      isMultiMediators <- TRUE
      num.meds <- length(fit.m)
    } else {
      isMultiMediators <- FALSE
      num.meds <- 1
    }
    ## Fit mediator model(s)
    # make an empty list to save mediator model(s)
    fit.meds <- vector(mode = "list", length = num.meds)
    for (i in 1:num.meds) {
      # The formula for each mediator model: lm(M1 ~ R + C); lm(M2 ~ R + C);...
      med.formula <- as.formula(paste(mediators[i], paste(c(treat, covariates), collapse = " + "), sep = " ~ "))
      if(is.factor(data[, mediators[i]])){
        fit.meds[[i]] <- lm(med.formula, data = data)
      } else {
        fit.meds[[i]] <- glm(med.formula, data = data)
      }
    }
    #fit.meds.new1 <- update(fit.meds[[1]], data = data.new)
    #vcov(fit.meds.new1)["racesex41", "racesex41"]
    
  } else if (!is.null(fit.m) & is.null(mediators)){
    
    # Possibility of multiple mediators
    if(length(class(fit.m)) == 1 && class(fit.m) == "list"){
      isMultiMediators <- TRUE
      num.meds <- length(fit.m)
    } else {
      isMultiMediators <- FALSE
      num.meds <- 1
    }
    
    fit.meds <- vector(mode = "list", length = num.meds)
    mediators <- rep(NA, num.meds)
    for (i in 1:num.meds) {
      fit.meds[[i]] <- fit.m[[i]]
      mediators[i] <- all.vars(formula(fit.m[[i]]))[[1]]
    }
    
  } else if (!is.null(fit.m) & !is.null(mediators)) {
    
    # Possibility of multiple mediators
    if(length(class(fit.m)) == 1 && class(fit.m) == "list"){
      isMultiMediators <- TRUE
      num.meds <- length(fit.m)
    } else {
      isMultiMediators <- FALSE
      num.meds <- 1
    }
    
    fit.meds <- vector(mode = "list", length = num.meds)
    mediators0 <- rep(NA, num.meds)
    for (i in 1:num.meds) {
      if(!isMultiMediators){
        fit.meds[[i]] <- fit.m
        mediators0[i] <- all.vars(formula(fit.m))[[1]]
      } else {
        fit.meds[[i]] <- fit.m[[i]]
        mediators0[i] <- all.vars(formula(fit.m[[i]]))[[1]]
      }
      if(mediators0[i] != mediators[i]){
        stop("Response variable(s) of 'fit.m' and 'mediators' must be match.")
      }
    }
    
  } else if (is.null(fit.m) & is.null(mediators)) {
    stop("Either 'fit.m' or 'mediators' must not be NULL.")
  }
  
  ##1. Calculate the SE of gamma_dm and df from the coefficient of M for a comparison group in fit.y.
  if (num.meds == 1) {
    meds <- all.vars(formula(fit.m))[[1]]
    if(is.factor(data[, mediators])){
      meds <- which(meds == substring(colnames(vcov(fit.y.new)), 1, nchar(meds)))[1]
    }
    var_gamma <- vcov(fit.y.new)[meds, meds]
    beta.resm.est <- coef(fit.y.new)[meds]
  } else if (num.meds == 2) {
    meds1 <- all.vars(formula(fit.m[[1]]))[[1]]
    meds2 <- all.vars(formula(fit.m[[2]]))[[1]]
    if(is.factor(data[, mediators[1]])){
      meds1 <- which(meds1 == substring(colnames(vcov(fit.y.new)), 1, nchar(meds1)))[1]
    }
    if(is.factor(data[, mediators[2]])){
      meds2 <- which(meds2 == substring(colnames(vcov(fit.y.new)), 1, nchar(meds2)))[1]
    }
    var_gamma <- vcov(fit.y.new)[meds1, meds1] + vcov(fit.y.new)[meds2, meds2] +
      2 * vcov(fit.y.new)[meds1, meds2]
    beta.resm.est <- sum(coef(fit.y.new)[c(meds1, meds2)]) ##??? when there are two Ms.
  } else if (num.meds > 2) {
    # will be added.
    stop("The case of three or more mediators is not supported yet.")
  }
  se_gamma <- sqrt(var_gamma)
  df <- fit.y.new$df.residual
  
  ##2. Calculate the effect of R on D and M (the coefficient for R in fit.m).
  treat.lev <- paste(treat, sel.lev.treat, sep = "")
  # sample covariance of initial disparity and disparity reduction
  cov.ini.red <- cov(boot.res$all.result[3 * (lev - 1) - 2, ],
                     boot.res$all.result[3 * (lev - 1), ])
  # sample covariance of initial disparity and alpha_r
  cov.ini.alphar <- cov(boot.res$all.result[3 * (lev - 1) - 2, ],
                        boot.res$alpha.r[, lev - 1])
  # sample covariance of initial disparity and se_gamma
  cov.ini.segamma <- cov(boot.res$all.result[3 * (lev - 1) - 2, ],
                         boot.res$se.gamma[, lev - 1])
  # sample variance of alpha_r
  var_alphahat.r <- var(boot.res$alpha.r[, lev - 1]) # 09/05/2022, 2/28/2023
  
  if (num.meds == 1) {
    alpha.r <- coef(fit.meds[[1]])[treat.lev] # coef of R[lev] on M = abs(alpha_r)
    ab <- abs(alpha.r)
    var_ab <- vcov(fit.meds[[1]])[treat.lev, treat.lev]
  } else if (num.meds == 2) {
    alpha.r <- coef(fit.meds[[1]])[treat.lev] + coef(fit.meds[[2]])[treat.lev]  # 09/05/2022
    ab <- abs(alpha.r)
    var_ab <- vcov(fit.meds[[1]])[treat.lev, treat.lev] +
      vcov(fit.meds[[2]])[treat.lev, treat.lev] - 2 * cov.ini.red
  } else if (num.meds > 2) {
    # will be added.
    stop("The case of three or more mediators is not supported yet.")
  }
  
  ## Calculate the variance of initial disparity # Use the one from the estimator function. 
  # Y ~ R + C
  ini.formula <- as.formula(paste(outcome, paste(c(treat, covariates), collapse = " + "), sep = " ~ "))
  ini <- lm(ini.formula, data = data)
  #var_tau <- vcov(ini)[treat.lev, treat.lev]
  var_tau <- var(boot.res$all.result[3 * (lev - 1) - 2, ])  # 09/05/2022
  
  ##3.
  # Point Estimate of Disparity Reduction for ref vs lev
  res.red <- boot.res$result[3 * (lev - 1), 1]
  # Point Estimate of Disparity Remaining for ref vs lev
  res.rem <- boot.res$result[3 * (lev - 1) - 1, 1]
  red <- rem <- lower_red <- lower_rem <- upper_red <- upper_rem <- var_ngamma <- matrix(NA, length(ry), length(rm))
  
  for (i in 1:length(ry)){
    for (j in 1:length(rm)){
      
      # true disp reduction (= est - bias)
      red[i, j] <- res.red - ab * se_gamma * sqrt(ry[i] * rm[j] * df)/sqrt(1 - rm[j]) #res.red
      beta.m <- beta.resm.est + ifelse(beta.resm.est > 0, - 1, 1) * se_gamma * sqrt(ry[i] * rm[j] * df)/sqrt(1 - rm[j]) # 09/05/2022
      # new variance of gamma
      var_ngamma[i, j] <- var_gamma * (1 - ry[i])/(1 - rm[j]) * (df/(df - 1)) 
      # new CI for disp reduction
      qn <- qnorm(1/2 + conf.level/2)
      se_red <- sqrt(alpha.r^2 * var_ngamma[i, j] + beta.m^2 * var_alphahat.r)  # Eq (15)
      lower_red[i, j] <- red[i, j] - qn * se_red
      upper_red[i, j] <- red[i, j] + qn * se_red
      
      # true disp remaining
      rem[i, j] <- res.rem - ab * se_gamma * sqrt(ry[i] * rm[j] * df)/sqrt(1 - rm[j]) #res.rem
      # new CI for disp remaining
      qn <- qnorm(1/2 + conf.level/2)
      k <- 1 * sqrt(ry[i] * rm[j] * df)/sqrt(1 - rm[j])
      se_rem <- sqrt(var_tau + se_red^2 - 2 * cov.ini.red +
                       2 * k * mean(boot.res$se.gamma[, lev - 1]) * cov.ini.alphar + # se_gamma
                       2 * k * mean(boot.res$alpha.r[, lev - 1]) * cov.ini.segamma)
      lower_rem[i, j] <- rem[i, j] - qn * se_rem
      upper_rem[i, j] <- rem[i, j] + qn * se_rem
      
    }
  }
  
  # output
  output <- c(lower_rem, upper_rem, lower_red, upper_red, red, rem)
  names(output) <- c("lower_rem", "upper_rem", "lower_red", "upper_red","red","rem")
  return(output)
  
}

# Coverage rate computation
cvgRate <- function(red.true, rem.true, red.mat, rem.mat){
  
  n.sim <- nrow(rem.mat)
  # coverage rate
  red.cvg <- rem.cvg <- 0
  for(i in 1:n.sim){
    red.cvg <- red.cvg + (red.true >= red.mat[i, 2] & red.true <= red.mat[i, 3])
    rem.cvg <- rem.cvg + (rem.true >= rem.mat[i, 2] & rem.true <= rem.mat[i, 3])
  }
  out <- matrix(c(red.cvg, rem.cvg), 1, 2)
  colnames(out) <- c("red.cvg", "rem.cvg")
  return(out/n.sim)
  
}
