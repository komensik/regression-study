##################################################
## POLISCI 672 Lab: Power and Covariates
## 2026-09-25
## Chern Xun Gan
##################################################


########################################
## Calculating Power
########################################

## STEP 1: GENERATE SAMPLE

# load required packages
library(tidyverse)
library(estimatr)
library(randomizr)

# set seed
set.seed(20260925)

# set number of obs
N <- 100

# binary indicator for whether an observation is a graduate (1) or undergraduate (0)
## 60% chance that a student is a graduate student
## each obs is a single draw from a binomial distribution with p = 0.6
matcha <- data.frame(grad = rbinom(N, 1, .6)) #drawing from binomial distribution

# for sample variability, don't want exactly 40 undergrads and 60 grads
# rbinom is random draw from binomial distribution
# 3 arguments: n = number observations, size = number of trials, prob is prob success on each trial 
# here is probability of being a graduate student (0.6)

# after doing this, get a row 



## STEP 2: SIMULATE TREATMENT ASSIGNMENT

# block random assignment
matcha <- matcha %>% 
  mutate(treat = block_ra(grad, .5))

# check treatment assignment
## there should be approximately 40% undergrads and 60% grads
## each group should have half of the subjects each in treatment and control

matcha %>% 
  group_by(grad) %>%
  summarize(n = n(), treated = sum(treat)) %>%
  ungroup()


## STEP 3: SIMULATE OUTCOME MEASURE
# different ways to do this step
# should be taking draws from binomial distribution

# generate outcomes
matcha <- matcha %>%
  mutate(
    # first assign probability of A for each subject based on DGP #
    prob_A = case_when( #probabilities are not the same for each row, so we give vector of probabilities
      grad == 0 & treat == 0 ~ .6, 
      grad == 0 & treat == 1 ~ .75, 
      grad == 1 & treat == 0 ~ .7, 
      grad == 1 & treat == 1 ~ .8
    ), 
    # then use probabilities to simulate exam score using random draw from binomial distribution
    exam = rbinom(N, 1, prob_A)
  )


## STEP 4: ESTIMATE ATE AND CHECK SUCCESS
# run model
model <- difference_in_means(exam ~ treat, matcha, grad)

# check if p-value < 0.05
success <- model$p.value < .05


## PUTTING IT ALL TOGETHER
# Set number of obs
N <- 100

for (i in 1:1000) {
  # binary indicator for whether an observation is a graduate (1) or undergraduate (0)
  matcha <- data.frame(grad = rbinom(N, 1, .6))
  
  # block random assignment
  matcha <- matcha %>% 
    mutate(treat = block_ra(grad, .5))
  
  # generate outcomes
  matcha <- matcha %>%
    mutate(
      # first assign probability of A for each subject based on DGP
      prob_A = case_when(
        grad == 0 & treat == 0 ~ .6, 
        grad == 0 & treat == 1 ~ .75, 
        grad == 1 & treat == 0 ~ .7, 
        grad == 1 & treat == 1 ~ .8
      ), 
      # then use probabilities to simulate exam score using random draw from binomial distribution
      exam = rbinom(N, 1, prob_A)
    )
  
  # run model
  model <- difference_in_means(exam ~ treat, matcha, grad)
  
  # check if p-value < 0.05
  success[i] <- model$p.value <= .05
}

# power for N = 100
mean(success)


## COMPARING POWER BETWEEN SAMPLE SIZES
powercalc <- function (N) {
  for (i in 1:1000) {
    # binary indicator for whether an observation is a graduate (1) or undergraduate (0)
    matcha <- data.frame(grad = rbinom(N, 1, .6))
    
    # block random assignment
    matcha <- matcha %>% 
      mutate(treat = block_ra(grad, .5))
    
    # generate outcomes
    matcha <- matcha %>%
      mutate(
        # first assign probability of A for each subject based on DGP
        prob_A = case_when(
          grad == 0 & treat == 0 ~ .6, 
          grad == 0 & treat == 1 ~ .75, 
          grad == 1 & treat == 0 ~ .7, 
          grad == 1 & treat == 1 ~ .8
        ), 
        # then use probabilities to simulate exam score using random draw from binomial distribution
        exam = rbinom(N, 1, prob_A)
      )
    
    # run model
    model <- difference_in_means(exam ~ treat, matcha, grad)
    
    # check if p-value < 0.05
    success[i] <- model$p.value <= .05
  }
  
  # power for sample size = N
  mean(success)
}

# N = 100, 150, 200, ..., 500
powers <- data.frame(sample_size = seq(100, 500, 50))

# calculate power for varying sample sizes
powers <- powers %>%
  mutate(power = sapply(sample_size, powercalc))

# plot power calculations
ggplot(powers, aes(sample_size, power)) + 
  geom_point() + 
  geom_smooth()

#############################
## USING DECLAREDESIGN
###########################

library(DeclareDesign)
#useful to know tidyverse to know how works, uses tidyverse

# declare initial sample size
sample_size = 100

# declare design
design <- 
  # declare the model (DGP), including sample size and the outcome
  declare_model(
    N = sample_size, # sample size
    X = rbinom(N, 1, .6), # covariate (block)
    potential_outcomes(Y ~ rbinom(N, 1, .6 + .1 * X + .15 * Z - .05 * X * Z)) # outcome
  ) + 
  # declare estimand (ATE in this case)
  declare_inquiry(ATE = mean(Y_Z_1 - Y_Z_0)) + 
  # declare random assignment (block random assignment in this case)
  declare_assignment(Z = block_ra(X, .5)) + 
  # declare measurement (telling DeclareDesign that the outcome can now be simulated)
  declare_measurement(Y = reveal_outcomes(Y ~ Z)) + #can p much just copy/paste this line onto any dec des that you do
  # declare estimator (block-adjusted DIM)
  declare_estimator(
    Y ~ Z, # give the equation
    inquiry = "ATE", 
    .method = difference_in_means, # specify which function (from estimatr) to use for estimator
    blocks = X # specify additional arguments for estimator function
  )

# use `redesign()` to change specific features of design
designs <- redesign(design, sample_size = seq(100, 500, 50))

# use `diagnose_design()` for power calculations (diagnosis test)
diagnosis <-
  diagnose_designs(
    designs, 
    diagnosands = declare_diagnosands(power = mean(p.value <= 0.05)), #power
    sims = 1000 # number of simulations for each sample size 
  )

# plot power calculation
ggplot(diagnosis$diagnosands_df, aes(x = sample_size, y = power)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = power - 1.96 * `se(power)`, ymax = power + 1.96 * `se(power)`)  
  ) + 
  geom_smooth(se = FALSE)


########################################
## Randomization checks
########################################
# import data
rush <- read.csv("Gendelman_2004.csv")

# Step 1: obtain test statistic (F-stat)
fstat <- lm_robust(treatment ~ pretest, rush)$fstatistic[1]

# Step 2: simulate random assignments
# declare random assignment
declaration <- declare_ra(nrow(rush), m = sum(rush$treatment))

# obtain permutation matrix
D <- obtain_permutation_matrix(declaration)

# initialize empty vector
f_sims <- numeric(ncol(D))

# randomization inference
for(i in 1:ncol(D)) {
  # Step 4: obtain "observed" data
  rush_ri <- rush %>% 
    mutate(d_temp = D[, i])
  
  # Step 5: obtain F-statistic from simulation
  f_sims[i] <- lm_robust(d_temp ~ pretest, rush_ri)$fstatistic[1]
}

# compute p-value
mean(f_sims >= fstat)

hist(f_sims, breaks = 100)
abline(v = fstat)

