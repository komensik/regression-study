library(readstata13)
library(boot)
library(marginaleffects)
library(data.table)


# Data from Cattaneo, Econometrica
dat <- read.dta13("birthwgtex.dta", convert.factors = F)
dat <- as.data.table(dat)


# OPTION: single model with interactions
jmod <- lm(bweight ~ mbsmoke * (prenatal1 + mmarried + mage + fbaby), data = dat)
summary(jmod)
print(avg_comparisons(jmod, variables = "mbsmoke", vcov = "HC1"))


# OPTION: model both outcomes
dat <- as.data.table(dat)
dat0 <- dat[mbsmoke == 0, ]
dat1 <- dat[mbsmoke == 1, ]

#                  | Y(0) - No Smoking | Y(1) - Smoking |
# -----------------|-------------------|----------------|
# T = 0 (No smoke) | Observed          | ?              |
# T = 1 (Smoke)    | ?                 | Observed       |
# -----------------|-------------------|----------------|

lm0 <- lm(bweight ~ prenatal1 + mmarried + mage + fbaby, data = dat0)
lm1 <- lm(bweight ~ prenatal1 + mmarried + mage + fbaby, data = dat1)

# predict missing (potential) oucome values
dat0[, bweight_0 := bweight]
dat1[, bweight_1 := bweight]
dat0[, bweight_1 := predict(lm1, newdata = dat0)]
dat1[, bweight_0 := predict(lm0, newdata = dat1)]

#                  | Y(0) - No Smoking | Y(1) - Smoking |
# -----------------|-------------------|----------------|
# T = 0 (No smoke) | dat0$bweight_0    | dat1$bweight_0 |
# T = 1 (Smoke)    | dat0$bweight_1    | dat1$bweight_1 |
# -----------------|-------------------|----------------|

dat_comb <- rbind(dat0,dat1)
dat_comb[, te := bweight_1 - bweight_0]

# ATE
mean(dat_comb$te)
# ATT
mean(dat_comb$te[dat_comb$mbsmoke==1])


#  We could get s.e. using the delta method. or use the bootstrap:

boot_ra <- function(dat, i){
    sdat <- dat[i,]
    dat0 <- sdat[mbsmoke == 0, ]
    dat1 <- sdat[mbsmoke == 1, ]
    lm0 <- lm(bweight ~ prenatal1 + mmarried + mage + fbaby, data = dat0)
    lm1 <- lm(bweight ~ prenatal1 + mmarried + mage + fbaby, data = dat1)
    # predict missing (potential) oucome values
    dat0[, bweight_0 := bweight]
    dat1[, bweight_1 := bweight]
    dat0[, bweight_1 := predict(lm1, newdata = dat0)]
    dat1[, bweight_0 := predict(lm0, newdata = dat1)]
    dat_comb <- rbind(dat0, dat1)
    dat_comb[, te := bweight_1 - bweight_0]
    est <- c(mean(dat_comb$te), mean(dat_comb$te[dat_comb$mbsmoke == 1]))
    return(est)
}

set.seed(67542)
samps <- boot(data = dat, statistic = boot_ra, R = 1999, 
    parallel = "multicore", ncpus = 16)
print(samps)
# get confidence bands: Confint(samps, type = "perc")



