
# Simulation to assess RI procedure for conducting multiple hypothesis tests
# See footnote 12, Chapter 9.

library(ri)

outeriter <- 100 # number of replications of overall procedure

reject <- rep(NA,outeriter)

for (iter in 1:outeriter) {

N <- 100  # sample size of hypothetical dataset

Z <- sample(c(rep(1,N/2),rep(0,N/2)))    # random assignment

X1 <- c(rep(1,N/4),rep(0,3*N/4))         # create covariates
X2 <- c(rep(0,3*N/4),rep(1,N/4))

Y <- runif(N) + 100*Z + 100*X1 #+ Z*X1

numiter <- 100   # iterations for purposes of calculating p-values

#perms <- genperms(Z,maxiter=numiter)


Yr <- Y - (mean(Y[Z==1])-mean(Y[Z==0]))*Z

storeP <- rep(NA,numiter)

for (i in 1:numiter) {
	Zri <- sample(Z)
	#perms[,i]
	storeP[i] <- min(summary(lm(Yr~Zri*X1+Zri*X2))$coefficients[5:6,4])  # store the minimum p-value
}

targetP <- sort(storeP)[round(numiter*.05)]


reject[iter] <- (min(summary(lm(Y~Z*X1+Z*X2))$coefficients[5:6,4]) <= targetP)

cat(iter,"")
}

summary(reject)

