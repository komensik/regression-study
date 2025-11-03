
library(readstata13)
library(data.table)
library(car)
library(lmtest)
library(sandwich)



################################################################# 
# Simple DiD example using Card and Krueger (1994) data
#################################################################

d <- as.data.table(read.dta13("CK1994/CK1994.dta"))
d

# CK definition of FTE:
d[, fte := empft + nmgrs + 0.5*emppt]

# Cell means by state * time
d[, mean(fte, na.rm = T), by = c("state", "time")]

# or, equivalently, using regression
m <- lm(fte ~ state * time, data = d)
b <- coef(m)
b["(Intercept)"]
b["(Intercept)"] + b["state"]
b["(Intercept)"] + b["time"]
b["(Intercept)"] + b["time"] + b["state"] + b["state:time"]

# To get standard errors: 
# 
# calculate using sqrt(k'Vk), where k is indicator vector of desired coefficient constraints,
# V is vcov(m) , e.g.  
# k = c(1,0,0,0); V = vcov(m); sqrt(t(k) %*% V %*% k)
# 
# or use linearHypothesis() from car
lh <- lht(m, c("(Intercept)=0", 
               "(Intercept)+state=0",
               "(Intercept)+time=0",
               "(Intercept)+state+time+state:time=0"), 
        vcov = hccm(m, "hc2"))
attr(lh, "value")
sqrt(diag(attr(lh, "vcov")))
# 
# or use the wrapper function lh_robust from the estimatr package 
# (which simply calls lm_robust() and linearHypothesis)
# lh_robust(fte ~ state * time, data = d, se_type = "HC2",
#     linear_hypothesis = c("(Intercept)=0", 
#                           "(Intercept)+state=0",
#                           "(Intercept)+time=0",
#                           "(Intercept)+state+time+state:time=0"
#                           ))



# DiD estimate
# Y = D + T + D*T
# can use the linear model from above!
# 

# First differences 
# For state 0 (delta_t)
b["time"]
# For state 1 (delta_t + D)
b["time"] + b["state:time"]

# Difference of difference
b["state:time"]
coeftest(m)["state:time",]
coeftest(m, vcov = vcovHC(m, "HC2"))["state:time",]




################################################################# 
# Repeated observations
#################################################################

load("panel.Rds")
d <- as.data.table(d)
# 200 cases
# 10 time periods, treatment occurs at t=5 for treat==1 group
# (note: timing does not vary)

m <- lm(y ~ treat + as.factor(id) + as.factor(time), data = d)
coeftest(m, vcov = vcovCL(m, cluster = ~ id))["treat",]


################################################################# 
# Test of parallelity of linear pre-trends
#################################################################

# Study trends in pre-treatment period
# Create variable indicating (time-invariant) treatment status
d[, tm := mean(treat), by = id]
d[tm == 0.5, C := 1]
d[tm == 0.0, C := 0]
pdf()
trends <- d[, mean(y), by = c("C", "time")]
with(trends, {
    plot(V1[C==0] ~ time[C==0], type = "l", col = "red", ylim = c(6, 8))
    lines(V1[C==1] ~ time[C==1], type = "l", col = "blue")
    })
abline(v = 5)
dev.off()


# identify pre- and post-treatment periods
d$dt0 <- ifelse(d$time<=5, 1, 0)
d$dt1 <- ifelse(d$time>5, 1, 0)

# Augmented model
pt <- lm(y ~ treat + as.factor(id) + as.factor(time) + time:dt0:C + time:dt1:C, data = d)
# Null hypothesis: parallel pre-trend
linearHypothesis(pt, "time:dt0:C=0")



#  Add leads to study anticipation effects

# create relative time to treatment 
d[, rtime := time - 5]

# manually create three pre-treatment lead indicators 
d[, lead1 := as.integer(rtime == -1 & C == 1)]
d[, lead2 := as.integer(rtime == -2 & C == 1)]
d[, lead3 := as.integer(rtime == -3 & C == 1)]

# more elegantly in one line...
# d[, paste0("lead", 1:3) := shift(y, n = 1:3, type = "lead"), by = id]

m2 <- lm(y ~ treat + lead3 + lead2 + lead1 + as.factor(id) + as.factor(time), data = d)
ct <- coeftest(m2, vcov = vcovCL(m2, cluster = ~ id))

# lead coefficients 
print(ct[c("lead3", "lead2", "lead1"), ])
# Joint test: all lead coefficients equal zero
linearHypothesis(m2, c("lead3=0", "lead2=0", "lead1=0"), vcov = vcovCL(m2, cluster = ~ id))



# Add lags and leads to study over-time changes
# sometimes called event-study 

# Load the event study example data
d <- readRDS("event_study_ex.rds")
d <- as.data.table(d)

# Define time grid (relative time periods to study)
event_grid <- -3:3

# Create indicator for treated units
d[, C := max(treated), by = id]

# Create event dummies for regression (ev_-3 .. e+v3) 
for(r in event_grid){
    nm <- if(r < 0) paste0("ev_m", abs(r)) else paste0("ev", r)
    d[, (nm) := as.integer(rtime == r & C == 1)]
}

# Include all event dummies and remove the intercept 
ev_names <- c("ev_m3", "ev_m2", "ev_m1", "ev0", "ev1", "ev2", "ev3")
f <- as.formula(paste0("ysim ~ 0 + ", paste(ev_names, collapse = " + "), " + as.factor(id) + as.factor(time)"))
mev <- lm(f, data = d)
ct_ev <- coeftest(mev, vcov = vcovCL(mev, cluster = ~ id))

# extract estimates and s.e. errors for the 7 periods
est <- rep(NA, length = 7)
se <- rep(NA, length = 7)
names(est) <- names(se) <- as.character(event_grid)
for (i in 1:length(event_grid)) {
    nm <- if (event_grid[i] < 0) {
            paste0("ev_m", abs(event_grid[i]))}
        else {
            paste0("ev", event_grid[i])
        }
    est[i] <- ct_ev[nm, "Estimate"]
    se[i] <- ct_ev[nm, "Std. Error"]
}


# PLOT
x <- 1:length(event_grid)
ylow <- est - 1.96 * se
yhigh <- est + 1.96 * se

par(mar = c(4,5,3,1))
plot(x, est, pch = 19, xlab = "Month", ylab = expression(beta[t]),
         ylim = c(min(ylow), max(yhigh)),
         cex.lab = 1.2, cex.axis = 1.1)
grid(nx = NA, ny = NULL, lty = 3, col = "lightgray")
abline(h = 0, col = "red")
arrows(x, ylow, x, yhigh, angle = 90, code = 3, length = 0.05)
points(x, est, pch = 19)


# Joint test for pre-trends (all leads = 0): 
lead_vars <- c("ev_m3", "ev_m2", "ev_m1")
linearHypothesis(mev, paste0(lead_vars, "=0"), vcov = vcovCL(mev, cluster = ~ id))


# use within estimator instead of adding i, t dummy variables
library(fixest)
f <- as.formula(paste0("ysim ~ 0 + ", paste(ev_names, collapse = " + "), " | id + time"))
mev <- fixest::feols(f, data = d, cluster = "id")
summary(mev)





#################################################################
# Using the DiD package
#################################################################

# Estimate group-time treatment effects

library(did)

d <- readRDS("did_ex.Rds")
head(d)
attgt <- att_gt(
    yname = "Y",
    tname = "period",
    idname = "id",
    gname = "G",
    data = d,
    est_method = "reg"
)
summary(attgt)
ggdid(attgt)

# conditional parallel trends
attgt <- att_gt(
    yname = "Y",
    tname = "period",
    idname = "id",
    gname = "G",
    data = d,
    est_method = "reg",
    xformla = ~X
)



# weighted average of all group-time effects
# but note, this over-weights effects of early treated groups
attgt_avg <- aggte(attgt, type = "simple")
summary(attgt_avg)

# Event-study style aggregation
attgt_es <- aggte(attgt, type = "dynamic")
summary(attgt_es)
ggdid(attgt_es)

# group-specific effects
attgt_gs <- aggte(attgt, type = "group")
summary(attgt_gs)
ggdid(attgt_gs)
