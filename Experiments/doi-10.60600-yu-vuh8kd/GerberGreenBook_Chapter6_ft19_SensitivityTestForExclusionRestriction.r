#-----------------------------------------------------------------------------
# Simulation to conduct a sensitivity test on violations of the exclusion restriction
# code by:       Al Fang, ahf2116@columbia.edu
# last updated:  26 jul 2013

# Code follows the model presented in
# Timothy G. Conley, Christian B. Hansen, and Peter E. Rossi.  2011.
# The Review of Economics and Statistics, February 2012, 94(1): 260–272.  
# See especially p.262.



#-----------------------------------------------------------------------------

rm(list=ls(all=TRUE))

library(AER)


#--------- define function to conduct simulations -----------------
# function name: ersim (exclusion restriction simulation)
#
# inputs:
#	g0 = vector of gamma_0 values in support G 
#			where no element is equal to zero
#	y = vector of outcome variable Y
#	x = vector of endogenous regressor
#	z = vector of instrument 
#
# outputs:
#	out = a matrix with cells containing b.hat, se.hat, p-values for each gamma0
#		from a 2SLS regression of (Y - Z*gamma0) on X using Z as instruments 

ersim <- function(g0,y,x,z){

# assign inputs to objects
Y <- as.vector(y)
X <- as.vector(x)
Z <- as.vector(z)
gamma0 <- as.vector(g0)

# create object to store results
out <- matrix(NA, nrow=length(gamma0), ncol=4) 
out[,1] <- gamma0
colnames(out) <- c("gamma0","b.hat","se.hat","p.value")

# loop over each value of gamma

for(i in 1:length(gamma0)){
	
	ystar <- Y-(gamma0[i] * Z)  # create new DV (Y - Z*gamma)
	fit <- ivreg(ystar ~ X | Z) # run 2SLS on the new DV
	b.hat <- summary(fit)$coefficients[2,1]   # pull estimated b
	se.hat <- summary(fit)$coefficients[2,2]  # pull estimated SE
	p.value <- summary(fit)$coefficients[2,4] # pull the p value
	out[i,2] <- b.hat
	out[i,3] <- se.hat
	out[i,4] <- p.value
	
}

return(out)

}


#-----------------------------------------------------------------------------
# EXAMPLE 1 with simulated data

# dgp

N <- 500
e1 <- rnorm(N, 0, 1)
e2 <- rnorm(N, 0, sd(e1)/3)
g <- rnorm(N, 5, 500)  		# coef on Z (2nd stage) -- true value, let gamma_0 vary in simulation under a constant effects assumption that gamma = gamma0
beta <-	rnorm(N, 3, 3)		# coef on X (2nd stage)
pi <- rnorm(N, 10, 5)		# coef on Z (1st stage)
z <- rnorm(N, 0, 1)
x <- .7 + z*pi + e1
y <- .5 + x*beta + z*g + e2	
data <- data.frame(y=y, x=x, z=z)
head(data)


# define set of gamma0's in support G

gamma.input <- seq(1,20,1)

# run function, returns p values for each gamma0

ersim(g0=gamma.input, y=data$y, x=data$x, z=data$z)


#-----------------------------------------------------------------------------
# EXAMPLE 2 with data from Miguel, Satyanath, Sergenti (2004)
#
# Full citation:
# Economic Shocks and Civil Conflict: An Instrumental Variables Approach
# Author(s): Edward Miguel, Shanker Satyanath, Ernest Sergenti
# Journal of Political Economy, Vol. 112, No. 4 (August 2004), pp. 725-753

setwd("/Users/donaldgreen/Dropbox/Field Experimentation Book/Datasets for Website/Sensitivity analysis -- chapter 7/")
library(foreign)
dat2 <- read.dta("mss_repdata.dta")
head(dat2)

# let y = any_prio_off
# let x = GPCP_g
# let z = gdp_g

# define set of gamma0's in support G

gamma.input <- c(seq(-2,-.01,.01), seq(.01, 2, .01))

# run function, returns p values for each gamma0

ersim(g0=gamma.input, y=dat2$any_prio_off, x=dat2$GPCP_g, z=dat2$gdp_g)
