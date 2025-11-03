
set.seed(67541)

library(readstata13)
library(MatchIt)
library(cobalt)
library(lmtest)
library(sandwich)



d <- read.dta13("lalonde.dta")
head(d)

# Some manual imbalance checking [canned versions below]
mean(d$hisp[d$treat == 1])
mean(d$hisp[d$treat == 0])

plot(density(d$educ[d$treat == 0], bw = 1))
lines(density(d$educ[d$treat == 1], bw = 1), col = "blue")


f <- treat ~ age + educ + black + hisp + married + nodegr


# Calculate propensity scores
ps_model <- glm(f, data = d, family = binomial(link = "probit"))
d$ps <- predict(ps_model, type = "response")

# and visually examine of overlap
plot(density(d$ps[d$treat == 0]), main = "Propensity Score Overlap", 
    xlab = "Propensity Score", col = "red", lwd = 2)
lines(density(d$ps[d$treat == 1]), col = "blue", lwd = 2)
legend("topright", c("Control", "Treatment"), col = c("red", "blue"), lwd = 2)


# Check for common support
min_treat <- min(d$ps[d$treat == 1])
max_treat <- max(d$ps[d$treat == 1])
min_control <- min(d$ps[d$treat == 0])
max_control <- max(d$ps[d$treat == 0])

cat("PS range treatment:", round(min_treat, 3), round(max_treat, 3), "\n")
cat("PS range control:  ", round(min_control, 3) , round(max_control, 3), "\n")



# MATCH. 1-NN using propensity score (probit of P(D|S))
m1 <- matchit(f,
    data = d, method = "nearest", distance = "glm", link = "probit")
# Study balance
summary(m1)
# PS balance
plot(m1, type = "jitter", interactive = FALSE)
bal.plot(m1)
# covariates
set.cobalt.options(binary = "std")
love.plot(m1, abs = T, drop.distance = TRUE, var.order = "unadjusted")

# MATCH. optimal matching on PS
m2 <- matchit(f,
    data = d, method = "optimal", distance = "glm", link = "probit")
summary(m2, un = F)
love.plot(m2, abs = T, drop.distance = TRUE, var.order = "unadjusted")
plot(m2, type = "jitter", interactive = FALSE)

# MATCH. Full matching 
m3 <- matchit(f, data = d, method = "full", distance = "glm", link = "probit")
summary(m3, un = F)
love.plot(m3, abs = T, drop.distance = TRUE, var.order = "unadjusted")
plot(m3, type = "jitter", interactive = FALSE)

# Create matched data sets
m1d <- match.data(m1)
m2d <- match.data(m2)
m3d <- match.data(m3)


# Estimate via OLS
# Adjust s.e. for pairs clustering

# M1 [Note: w_ij=1 in this case (1-NN match)]
m1_ols <- lm(re78 ~ treat, data = m1d, weights = weights)
coeftest(m1_ols, vcov. = vcovCL, cluster = ~subclass)

# M3
m3_ols <- lm(re78 ~ treat, data = m3d, weights = weights)
coeftest(m3_ols, vcov. = vcovCL, cluster = ~subclass)

# Include covariates in outcome regression [remaining imbalances]
m3_ols_cov <- lm(re78 ~ treat + age + educ + black + hisp + 
    married + nodegr, 
    data = m3d, weights = weights)
coeftest(m3_ols_cov, vcov. = vcovCL, cluster = ~subclass)




# Use DISTANCE measures other than PS
# Full matching using MH distance
m4 <- matchit(f, data = d, method = "full", distance = "mahalanobis")
summary(m4, un = F)
m4d <- match.data(m4)
m4_ols <- lm(re78 ~ treat, data = m4d, weights = weights)
coeftest(m4_ols, vcov. = vcovCL, cluster = ~subclass)



