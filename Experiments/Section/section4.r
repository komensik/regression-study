##################################################
## POLISCI 672 Lab: Randomization Inference
## 2026-09-18
## Chern Xun Gan
##################################################


########################################
## Blocking and clustering
########################################

## DECLARING RANDOM ASSIGNMENT

declaration <- 
  declare_ra(
    nrow(camerer), 
    camerer$pair, 
    block_m = rep(1, length(unique(camerer$pair)))
  )

declaration1 <- 
  declare_ra(nrow(camerer), clusters = camerer$pair, m = 8)

declaration2 <- 
  declare_ra(nrow(camerer), camerer$venue, camerer$pair, block_m = c(4, 4))



camerer <- 
  camerer %>% 
  mutate(d = conduct_ra(declaration), 
         d1 = conduct_ra(declaration1), 
         d2 = conduct_ra(declaration2))


## ESTIMATING ATE

# block
camerer %>% 
  difference_in_means(experimentbets ~ d, ., pair)

# cluster
camerer %>% 
  difference_in_means(experimentbets ~ d1, ., cluster = pair)

# block and cluster
camerer %>% 
  difference_in_means(experimentbets ~ d2, ., venue, pair)



# block
camerer %>% 
  lm_robust(experimentbets ~ d, ., fixed_effects = pair)

# cluster
camerer %>% 
  lm_robust(experimentbets ~ d1, ., clusters = pair)

# block and cluster
camerer %>% 
  lm_robust(experimentbets ~ d2, ., clusters = pair, fixed_effects = venue)



camerer <- 
  camerer %>% 
  mutate(ipw = 1 / obtain_condition_probabilities(declaration2, d2))



camerer %>% 
  difference_in_means(experimentbets ~ d2, ., venue, pair)

camerer %>% 
  lm_robust(experimentbets ~ d2 + as.factor(venue), 
            ., 
            ipw, 
            clusters = pair
  )


########################################
## Confidence intervals from randomization inference
########################################

# Step 1: obtain ATE estimate
dim <- camerer %>%
  difference_in_means(experimentbets ~ treatment, ., pair) %>% 
  coefficients()

# Step 2: construct ATE estimates
camerer_ri <- 
  camerer %>% 
  mutate(
    y0 = case_when(
      treatment == 0 ~ experimentbets, 
      treatment == 1 ~ experimentbets - dim
    ), 
    y1 = case_when(
      treatment == 1 ~ experimentbets, 
      treatment == 0 ~ experimentbets + dim
    )
  )

# Step 3: simulate random assignments
D <- obtain_permutation_matrix(declaration)

# initialize empty vector
dim_sims <- c()

# obtain 10,000 DIM estimates from simulation
for (i in 1:10000) {
  # Step 4: obtain "observed" outcomes
  camerer_ri <- 
    camerer_ri %>% 
    mutate(
      d_temp = D[, i], 
      y_temp = case_when(
        d_temp == 0 ~ y0, 
        d_temp == 1 ~ y1
      )
    )
  
  # Step 5: compute DIM
  dim_sims[i] <- 
    camerer_ri %>% 
    difference_in_means(y_temp ~ d_temp, ., pair) %>% 
    coefficients()
}


# 95% confidence interval
camerer_ci <- dim_sims %>% 
  quantile(c(.025, .975))

hist(dim_sims, breaks = 100)
abline(v = dim, col = "blue")
abline(v = camerer_ci, col = "red")