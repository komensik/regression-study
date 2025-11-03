
library(readstata13)
library(MSCMT)
library(ggplot2)

set.seed(765222)


#####################################################################
# Read data
#####################################################################
d <- read.dta13("CaliProp99.dta") 

# Convert to list of TxN matrices
dl <- listFromLong(d, 
    unit.variable = "state_num",
    time.variable =  "year", 
    unit.names.variable = "state")



#####################################################################
# Basic SCM
#####################################################################

treated <- "California"
states <- colnames(dl[[1]])
controls  <- setdiff(states, treated)

# 2xL matrix of LHS variables, Y
times_dep <- cbind("cigsale" = c(1970, 1988))

# 2xK matrix of predictor variables, X
times_pred <- cbind(
    "age15to24" = c(1970, 1988),
    "retprice" = c(1970, 1988),
    "lnincome" = c(1970, 1988),
    "beer" = c(1984, 1988),
    "cigsale" = c(1975, 1988)
)

agg_fns <- c(rep("mean", 4), "id")

m1 <- mscmt(dl, treatment.identifier = treated, 
    controls.identifier = controls,
    times.dep = times_dep, times.pred = times_pred, 
    agg.fns = agg_fns)
print(m1)
plot(m1, type = "comparison")
plot(m1, type = "gaps")



#####################################################################
# Placebo treatments
#####################################################################

m1p <- mscmt(dl, treatment.identifier = treated, 
    controls.identifier = controls,
    times.dep = times_dep, times.pred = times_pred, 
    agg.fns = agg_fns, placebo = T)

par(mfrow = c(1,3))
plot(m1p[["Nevada"]], type = "comp", ylim = c(0, 300))
plot(m1p[["Kentucky"]], type = "comp", ylim = c(0, 300))
plot(m1, type = "comp", ylim = c(0, 300))

# Plot gaps and placebos
plot(m1p[[treated]]$gaps$cigsale, ylab = "Gap")
for(i in controls){
    lines(m1p[[i]]$gaps$cigsale, col = "gray60")
}

# Exclude cases w/ high rmspe compared to treated case
plot(m1p[[treated]]$gaps$cigsale, ylab = "Gap", lwd = 4)
for(i in controls){
    if(m1p[[i]]$rmspe < m1p[[treated]]$rmspe * 5 ) {
        lines(m1p[[i]]$gaps$cigsale, col = "gray60")
    } 
}
lines(m1p[[treated]]$gaps$cigsale, col = "red", lwd =2)
abline(v = 1988, lty = 2)


# easy-to-use version using ggplot:
# plot(m1p)
ggplot(m1p, exclude.ratio = 5, ratio.type="rmspe")


# Based on placebos:
# Per-period p-values
pvalue(m1p, exclude.ratio = 5, ratio.type="rmspe")
# DiD type overall pre-post difference
did(m1p, exclude.ratio = 5)




#####################################################################
# Good practice
# Multiple starting sets (via vector of RNG seeds)
#####################################################################
seeds <- sample(5555:87653, 20, replace = FALSE) 

# Optional: execute in parallel
library(parallel)
cl <- makeCluster(4)  #use whatever number of cores you'd like to use

m2 <- mscmt(dl, treatment.identifier = treated, 
    controls.identifier = controls,
    times.dep = times_dep, times.pred = times_pred, 
    agg.fns = agg_fns,
    seed = seeds, check.global = TRUE, placebo = T,
    cl = cl
)
print(m2)
plot(m2)




#####################################################################
# Placebo treatment 1988 ==> 1980
#####################################################################

times_dep <- cbind("cigsale" = c(1970, 1980))
times_pred <- cbind(
    "age15to24" = c(1970, 1980),
    "retprice" = c(1970, 1980),
    "lnincome" = c(1970, 1980),
    "cigsale" = c(1975, 1980)
)
agg_fns <- c(rep("mean", 3), "id")
m3 <- mscmt(dl, treatment.identifier = treated, 
    controls.identifier = controls,
    times.dep = times_dep, times.pred = times_pred, 
    agg.fns = agg_fns)
plot(m3, type = "comparison", legend = F)
abline(v=1980, lty = 3)
