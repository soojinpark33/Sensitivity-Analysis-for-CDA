# Load package
library(causal.decomp)

# Load data: Try '?sMIDUS' to get more information on the sMIDUS data.
data(sMIDUS)

# 'sMIDUS' is a synthetic dataset that includes variables from the Midlife Development
# in the U.S. (MIDUS) study, which is described in the manuscript. It has been artificially
# generated based on the actual MIDUS data, which is not publicly available due to
# confidentiality concerns.
#
# The tables in the manuscript cannot be replicated exactly using this synthetic data
# since it is artificially generated.

# Center covariates
sMIDUS$age <- scale(sMIDUS$age, center = TRUE, scale = FALSE)
sMIDUS$stroke <- scale(sMIDUS$stroke, center = TRUE, scale = FALSE)
sMIDUS$T2DM <- scale(sMIDUS$T2DM, center = TRUE, scale = FALSE)
sMIDUS$heart <- scale(sMIDUS$heart, center = TRUE, scale = FALSE)

#-----------------------------------------------------------------------------------------------#
# Table 1
#-----------------------------------------------------------------------------------------------#
# M ~ R + C
fit.m <- lm(edu ~ racesex + age + stroke + T2DM + heart, data = sMIDUS)
# Y ~ R + X + M + C
fit.y1 <- lm(health ~ racesex + lowchildSES + abuse + edu + age + stroke + T2DM + heart, data = sMIDUS)
# Y ~ R + M + R:M + C
fit.y2 <- lm(health ~ racesex + edu + racesex:edu + age + stroke + T2DM + heart, data = sMIDUS)
# Y ~ R + X + M + R:M + C
fit.y3 <- lm(health ~ racesex + lowchildSES + abuse + edu + racesex:edu + age + stroke + T2DM + heart, data = sMIDUS)

# First column: R-M interaction (N), intermediate confounders (Y)
tab1.1 <- smi(fit.m = fit.m, fit.y = fit.y1, sims = 1000, conf.level = .95,
              covariates = c("age", "stroke", "T2DM", "heart"), treat = "racesex", seed = 227)
tab1.1

# Second column: R-M interaction (Y), intermediate confounders (N)
tab1.2 <- smi(fit.m = fit.m, fit.y = fit.y2, sims = 1000, conf.level = .95,
              covariates = c("age", "stroke", "T2DM", "heart"), treat = "racesex", seed = 227)
tab1.2

# Third column: R-M interaction (Y), intermediate confounders (Y)
tab1.3 <- smi(fit.m = fit.m, fit.y = fit.y3, sims = 1000, conf.level = .95,
            covariates = c("age", "stroke", "T2DM", "heart"), treat = "racesex", seed = 227)
tab1.3

#-----------------------------------------------------------------------------------------------#
# Figure 2. Sensitivity analysis Using Regression Coefficients
#-----------------------------------------------------------------------------------------------#
len.treat <- length(levels(sMIDUS[, "racesex"]))
sel.lev.treat = "4" # To compare Black Women (R = 4) with White Men (R = 1)
lev <- which(levels(sMIDUS[, "racesex"]) == sel.lev.treat)
# Note: In the paper, Black Women (R = 1), White Men (R = 0)

# Point Estimate of Disparity Reduction for ref vs lev
res.red <- tab1.3$result[3 * (lev - 1), 1]
res.red.lo <- tab1.3$result[3 * (lev - 1), 2]
res.red.up <- tab1.3$result[3 * (lev - 1), 3]
# Point Estimate of Disparity Remaining for ref vs lev
res.rem <- tab1.3$result[3 * (lev - 1) - 1, 1]
res.rem.lo <- tab1.3$result[3 * (lev - 1) - 1, 2]
res.rem.up <- tab1.3$result[3 * (lev - 1) - 1, 3]

alpha_r <- coef(fit.m)[4]

# Create sensitivity contour plots
par(mfrow = c(1, 2))
eta.by = 0.05
delta <- seq(0 + eta.by, 1.5 - eta.by, 0.01)
beta  <- seq(0 + eta.by, 1.5 - eta.by, 0.01)
eff <- lower <- upper <- matrix(NA, length(delta), length(beta))

# Disparity reduction
for(i in 1:length(delta)){
  for(j in 1:length(beta)){
    eff[i, j] <- abs(res.red) - abs(alpha_r) * delta[i] * beta[j]
    lower[i, j] <- res.red.lo - alpha_r * delta[i] * beta[j]
    upper[i, j] <- res.red.up - alpha_r * delta[i] * beta[j]
  }
}
contour(beta, delta, eff, levels = 0, lwd = 2, main = "A) Disparity reduction",
        xlab = expression(beta[u]), ylab = expression(delta[m]), ylim = c(0, 1.5))
contour(beta, delta, eff, levels = seq(-3, 3, 0.2), lwd = 1, add = TRUE, lty = 1)

# Disparity remaining
for(i in 1:length(delta)){
  for(j in 1:length(beta)){
    eff[i, j] <- abs(res.rem) - abs(alpha_r) * delta[i] * beta[j]
    lower[i, j] <- res.rem.lo - alpha_r * delta[i] * beta[j]
    upper[i, j] <- res.rem.up - alpha_r * delta[i] * beta[j]
  }
}
contour(beta, delta, eff, levels = 0, lwd = 2, main = "B) Disparity remaining",
        xlab = expression(beta[u]), ylab = expression(delta[m]), ylim = c(0, 1.5))
contour(beta, delta, eff, levels = seq(-3, 3, 0.2), lwd = 1, add = TRUE, lty = 1)

#-----------------------------------------------------------------------------------------------#
# Figure 3. Sensitivity analysis Using R-Squared Values
#-----------------------------------------------------------------------------------------------#
# Note: The argument 'sel.lev.treat' specifies the level of treatment,
#       which we want to compare with the reference level.
# In our case, the reference level is "1" (White men) and we compare it to "4" (Black women). 
sensRes <- sensitivity(boot.res = tab1.3, fit.m = fit.m, fit.y = fit.y3, mediators = "edu",
                       covariates = c("age", "stroke", "T2DM", "heart"), treat = "racesex",
                       sel.lev.treat = "4", max.rsq = 0.3)
names(sensRes)
sensRes
sensRes$RV_red       # 0.1230061
sensRes$RV_red_alpha # computational estimate = 0.1212143
sensRes$RV_rem       # 0.08611574
sensRes$RV_rem_alpha # computational estimate = 0.0847163

# Create sensitivity contour plots
plot(sensRes)
