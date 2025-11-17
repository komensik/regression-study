
library(ivreg)
library(sandwich)
library(lmtest)

set.seed(836365)


#  ============================
#  =       Simulation         =
#  ============================

# Func: generate data following DAG used in the lecture
# Y, D, given Z, U
mkdata <- function(d.to.y, z.to.d, u.to.y, u.to.d, z.to.y, N){
  Z <- rnorm(N) ## instrument
  U <- rnorm(N) ## unobs. confounder 
  D <- Z * z.to.d + U * u.to.d + rnorm(N, sd=1)
  Y <- D * d.to.y + U * u.to.y + Z * z.to.y + rnorm(N, sd=1)
  data.frame(D = D, Z = Z, U = U, Y = Y)
}


# Data under valid exclusion, strong instr.
dat <- mkdata(d.to.y = 1, z.to.d = 0.5, u.to.y = 0.5, 
    u.to.d = 0.5, z.to.y = 0, N = 1000) 

# First stage
summary(lm(D ~ Z , data = dat))

# IV estimate
m1 <- ivreg(Y ~ D  | Z, data = dat)
summary(m1, diagnostics = TRUE)
coeftest(m1, vcov = vcovHC(m1, type = "HC1"))


# Quick Simulation
nsim <- 1000
sims <- rep(NA, nsim)
for(i in 1:nsim){
  d_sim <- mkdata(d.to.y = 1, z.to.d = 0.4, u.to.y = 0.5, 
      u.to.d = 0.5, z.to.y = 0, N = 1000) 
  m <- ivreg(Y ~ D | Z, data = d_sim)
  m <- coeftest(m, vcov = vcovHC(m, type = "HC1"))
  sims[i] <- m["D",1]
}
mean(sims) - 1


# weak instrument
dat <- mkdata(d.to.y = 1, z.to.d = 0.04, u.to.y = 0.5, 
    u.to.d = 0.5, z.to.y = 0, N = 1000) 

m2 <- ivreg(Y ~ D  | Z, data = dat)
summary(m2, diagnostics = TRUE)
coeftest(m2, vcov = vcovHC(m2, type = "HC1"))


# Simulation
sims <- rep(NA, nsim)
for(i in 1:nsim){
  d_sim <- mkdata(d.to.y = 1, z.to.d = 0.04, u.to.y = 0.5, 
      u.to.d = 0.5, z.to.y = 0, N = 1000) 
  m <- ivreg(Y ~ D | Z, data = d_sim)
  m <- coeftest(m, vcov = vcovHC(m, type = "HC1"))
  sims[i] <- m["D",1]
}
mean(sims) - 1


# violated exclusion 
dat <- mkdata(d.to.y = 1, z.to.d = 0.1, u.to.y = 0.5, 
    u.to.d = 0.5, z.to.y = 0.2, N = 1000) 

m3 <- ivreg(Y ~ D  | Z, data = dat)
summary(m3, diagnostics = TRUE)
coeftest(m3, vcov = vcovHC(m3, type = "HC1"))

# Simulation
sims <- rep(NA, nsim)
for(i in 1:nsim){
  d_sim <- mkdata(d.to.y = 1, z.to.d = 0.1, u.to.y = 0.5, 
      u.to.d = 0.5, z.to.y = 0.2, N = 1000) 
  m <- ivreg(Y ~ D | Z, data = d_sim)
  m <- coeftest(m, vcov = vcovHC(m, type = "HC1"))
  sims[i] <- m["D",1]
}
mean(sims) - 1


