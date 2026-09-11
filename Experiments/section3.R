##################################################
## POLISCI 672 Lab: Randomization Inference
## 2026-09-11
## Chern Xun Gan
##################################################

########################################
## A quick refresher of for-loops
########################################

# calculates values but does not print or store them
for (i in 1:3) {
  2 * i
}

# an output will be produced only if explicitly called
for (i in 1:3) {
  print(2 * i)
}

# however, we still cannot save the for-loop output to an object (check!)
loop_out <- 
  for (i in 1:3) {
    print(2 * i)
  }

# the simplest way to save your for-loop output is to first define some empty object
loop_out <- c()
# then save the output from each iteration into the object
for (i in 1:3) {
  loop_out[i] <- 2 * i
}

loop_out

# I sometimes like to create the empty object inside the for-loop
# this is to prevent adding to already-created objects
for (i in 1:3) {
  # check if empty vector should be created
  # note however that this forces the loop to make this check on every iteration
  if (i == 1) {
    loop_out <- c()
  }
  loop_out[i] <- 2 * i
}

loop_out


# initialize empty vector
loop_out <- c()

# the method of using `i == 1` to check in the loop will not work here
# this is because the first element of the input vector is not 1
for (i in 2:4) {
  loop_out[i] <- 2 * i
}

# the first element of `loop_out` is empty
loop_out

# `loop_out` will be overwritten with every iteration
for (i in 1:3) {
  loop_out <- 2 * i
}

loop_out


# `loop_out` already exists, and we are now adding to the existing vector
for (i in 1:3) {
  loop_out <- c(loop_out, 2 * i)
}

loop_out


########################################
## Random assignment with `randomizr`
########################################

# load required packages
library(tidyverse)
library(estimatr)
library(randomizr)

# import data
hajj <- read.csv("Clingingsmith_et_al_QJE_2009dta.csv")

# always set seed when using randomizr!
set.seed(20260911)


hajj_simple <- 
  simple_ra(nrow(hajj), .5)

hajj_complete_m <- 
  complete_ra(nrow(hajj), sum(hajj$success))

hajj_complete_p <- 
  complete_ra(nrow(hajj), prob = .5)

hajj_three_arm <- 
  complete_ra(nrow(hajj), prob_each = rep(1/3, 3), conditions = 0:2)


# `hajj_complete_m` should sum to the same as `success`
sum(hajj_complete_m) == sum(hajj$success)

# `hajj_complete_p` should have a mean of 0.5
mean(hajj_complete_p)

# `hajj_three_arm` should show a roughly equal distribution across 0, 1, 2
table(hajj_complete_m)



declaration <- declare_ra(nrow(hajj), prob = .5, simple = TRUE)

declaration


declaration %>% 
  conduct_ra() %>% 
  head()



D <- obtain_permutation_matrix(declaration)

# preview permutation matrix
D[1:5, 1:5]

# dimensions of permutation matrix
dim(D)


########################################
## Randomization inference with `randomizr` and for-loops
########################################

## STEP ONE
# using switching equation
hajj_ri <- 
  hajj %>% 
  mutate(
    y0_null = (1 - success) * views + success * (views - 0.25), 
    y1_null = success * views + (1 - success) * (views + 0.25)
  )

# using case_when
hajj_ri <- 
  hajj %>% 
  mutate(
    y0_null = case_when(
      success == 0 ~ views, 
      success == 1 ~ views - 0.25
    ), 
    y1_null = case_when(
      success == 1 ~ views, 
      success == 0 ~ views + 0.25
    )
  )

# preview data
hajj_ri %>% 
  select(success, views, y0_null, y1_null) %>%
  head()

# check if `y1_null - y0_null == .25` for each subject, then see if `TRUE` uniquely exists
unique(hajj_ri$y1_null - hajj_ri$y0_null == .25)


## STEP TWO
# create random assignment vector
hajj_ri <- 
  hajj_ri %>% 
  mutate(d_temp = complete_ra(nrow(hajj_ri), sum(success)))

# check if random assignment is correct
sum(hajj_ri$success) == sum(hajj_ri$d_temp)

# preview data
hajj_ri %>% 
  select(success, views, y0_null, y1_null, d_temp) %>%
  head()


## STEP THREE
# using switching equation
hajj_ri <- 
  hajj_ri %>% 
  mutate(y_temp = d_temp * y1_null + (1 - d_temp) * y0_null)

# using case_when
hajj_ri <- 
  hajj_ri %>% 
  mutate(
    y_temp = case_when(
      d_temp == 0 ~ y0_null, 
      d_temp == 1 ~ y1_null
    )
  )

# preview data
hajj_ri %>% 
  select(success, views, y0_null, y1_null, d_temp, y_temp) %>%
  head()


## STEP FOUR
# using DIM function
sharp_null_dim <- 
  hajj_ri %>% 
  difference_in_means(y_temp ~ d_temp, .) %>% 
  coefficients

# using lm_robust
sharp_null_dim <- 
  hajj_ri %>% 
  lm_robust(y_temp ~ d_temp, .) %>%
  coefficients %>%
  .[2]


## PUTTING IT ALL TOGETHER
sharp_null_dim <- c()

hajj_ri <- 
  hajj %>% 
  mutate(
    y0_null = case_when(
      success == 0 ~ views, 
      success == 1 ~ views - .25
    ), 
    y1_null = case_when(
      success == 1 ~ views, 
      success == 0 ~ views + .25
    )
  )

for (i in 1:10000) {
  hajj_ri <- 
    hajj_ri %>% 
    mutate(
      d_temp = complete_ra(nrow(hajj_ri), sum(success)), 
      y_temp = case_when(
        d_temp == 0 ~ y0_null, 
        d_temp == 1 ~ y1_null
      )
    )
  
  sharp_null_dim[i] <- 
    hajj_ri %>%
    lm_robust(y_temp ~ d_temp, .) %>% 
    coefficients %>%
    .[2]
}


length(sharp_null_dim)
hist(sharp_null_dim)
mean(sharp_null_dim)


declaration <- declare_ra(nrow(hajj), m = sum(hajj$success))

D <- obtain_permutation_matrix(declaration)


## CALCULATING P-VALUES
# compute DIM from actual experiment
hajj_dim <- 
  hajj %>% 
  difference_in_means(views ~ success, .) %>% 
  coefficients

hajj_dim

# number of simulated DIM at least as large as actual
sum(sharp_null_dim >= hajj_dim)

# one-sided p-value
mean(sharp_null_dim >= hajj_dim)

# number of simulated DIM at least as large in absolute value as actual
sum(abs(sharp_null_dim) >= abs(hajj_dim))

# two-sided p-value
mean(abs(sharp_null_dim) >= abs(hajj_dim))


hist(sharp_null_dim, xlim = c(- .5, .9))
abline(v = hajj_dim, col = "blue")
abline(v = - hajj_dim, col = "red")
