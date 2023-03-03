# Load package and source file
library(causal.decomp)
source("SensAnal_source.R")

#-----------------------------------------------------------------------------------------------#
# Appendix E: Simulations (Table 2)
#-----------------------------------------------------------------------------------------------#
# We are going to generate population data with different effect sizes, using sMIDUS data.
# First, pick two race-gender groups and change labels as given in the manuscript (just for consistency)
data(sMIDUS)
sMIDUS <- sMIDUS[sMIDUS$racesex %in% c("1", "4"), ] # 1: White men, 4: Black women
sMIDUS$racesex <- factor(sMIDUS$racesex)
levels(sMIDUS$racesex) <- c("0", "1") # 0: White men, 1: Black women

# Fit models on M and Y to extract regression coefficients from sMIDUS data.
age.grp <- cut(sMIDUS$age, breaks = c(0, 50, Inf), include.lowest = TRUE, labels = c(1, 2))
fit.X <- lm(lowchildSES ~ racesex + age.grp, data = sMIDUS)
fit.U <- lm(T2DM ~ racesex + lowchildSES + age.grp, data = sMIDUS)
fit.M <- lm(edu ~ racesex + lowchildSES + T2DM + age.grp, data = sMIDUS)
fit.Y <- lm(health ~ racesex + lowchildSES + edu + T2DM + age.grp, data = sMIDUS)
coef.X <- coef(fit.X) # X ~ R + C.grp
coef.U <- coef(fit.U) # U ~ R + X + C.grp
coef.M <- coef(fit.M) # M ~ R + X + U + C.grp
coef.Y <- coef(fit.Y) # Y ~ R + X + M + U + C.grp

# Generate population data with different effect sizes
pop.small <- popDataGen(effectSize = "small", seed = 10)    # so that ry = rm = 0.02
pop.medium <- popDataGen(effectSize = "medium", seed = 10)  # so that ry = rm = 0.13
pop.large <- popDataGen(effectSize = "large", seed = 10)    # so that ry = rm = 0.26

# In Appendix E, we consider nine different conditions;
# sample size: 100, 500, 1000; and effect size: small, medium, large.
# For example, try the following setting of n = 100 and effect sizes = small.
# The setting can be changed when needed.
n <- 100
pop.dat <- pop.small$data
ry <- pop.small$ry             # in this case, effect size = small, so ry = 0.02
rm <- pop.small$rm             # in this case, effect size = small, so rm = 0.02
Red.true <- pop.small$Red.true # true value of disparity reduction
Rem.true <- pop.small$Rem.true # true value of disparity remaining

n.sim <- 100 # number of simulations
n.cores <- 1  # number of cores used in 'smi' function.

rem.CI <- red.CI <- matrix(NA, n.sim, 3)
colnames(rem.CI) <- c("estimate", "lower", "upper")
colnames(red.CI) <- c("estimate", "lower", "upper")
for(i in 1:n.sim){
  
  set.seed(i)
  samp.ind <- sample(nrow(pop.dat), n, replace = FALSE)
  dat <- data.frame(pop.dat[samp.ind, ])
  # Center covariate(s)
  dat$C.grp <- relevel(dat$C.grp, ref = "2")
  
  fit.m <- lm(M ~ R  + C.grp , data = dat)
  fit.y <- lm(Y ~ R + X + M + C.grp, data = dat)
  smi.res <- smi(fit.m = fit.m, fit.y = fit.y, sims = 1000, conf.level = .95,
                 covariates = "C.grp", treat = "R", seed = i, mc.cores = n.cores)
  sens.res <- sens.for.se(boot.res = smi.res, fit.y = fit.y, fit.m = fit.m, mediators = "M",
                          covariates = "C.grp", treat = "R", sel.lev.treat = "1", ry = ry, rm = rm)
  
  rem.CI[i, ] <- c(smi.res$result[2, 1], sens.res[c(1, 2)]) # estimate, lower, upper
  red.CI[i, ] <- c(smi.res$result[3, 1], sens.res[c(3, 4)]) # estimate, lower, upper
  print(i)
  
}

# Compute coverage rates (Table 2 in Appendix E)
cvgRate(red.true = Red.true, # true value of disparity reduction
        rem.true = Rem.true, # true value of disparity remaining
        red.mat = red.CI,    # 95% CIs of disparity reduction, which are constructed using our SE formula
        rem.mat = rem.CI)    # 95% CIs of disparity remaining, which are constructed using our SE formula
