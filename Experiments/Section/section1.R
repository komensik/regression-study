##################################################
## POLISCI 672 Lab: Starting Up
## 2026-08-28
## Chern Xun Gan
##################################################

########################################
## Data transformation with `dplyr`
########################################

# load packages
library(tidyverse)

# load data
load("caffeine.Rdata")

# show first few rows
head(caffeine)


### Pipes ###

# this is what we would do in base R
output1 <- 1 + 1
output2 <- output1 + 1
rep(output2, 5)

# this is the equivalent when using pipes
1 + 1 %>%   ## the output for this line is 2, we pass this output to the next line
  + 1 %>%   ## therefore, here we are doing (1 + 1) + 1 = 3, and pass it on
  rep(5)    ## this is now equivalent to rep((1 + 1) + 1, 5)

2 %>% c(1) # with the pipe, "1" becomes the second element of the vector!


### Subsetting and creating new variables ###

# coerce to numeric and remove rows with missing data
caffeine <- 
  caffeine %>% 
  sapply(as.numeric) %>%
  as.data.frame() %>%
  na.omit()

# show observations that received the matcha treatment
caffeine %>% 
  filter(matcha == 1) %>%
  head()

# keep only coffee treatment status and score
caffeine %>%
  select(coffee, score) %>%
  head()

# remove observations that received the coffee treatment, remove the coffee 
# treatment variable, and assign this to a new object
nocoffee <- 
  caffeine %>% 
  filter(coffee != 1) %>% 
  select(!coffee)
## note that doing this assigns as an object the end result after all pipes!

# define a new variable `treated` as assigned to any treatment condition
caffeine <- 
  caffeine %>% 
  mutate(
    treat = 
      (matcha + coffee > 0) %>% 
      as.numeric() ## you can also pipe inside functions!
  )

head(caffeine)
## it is always a good idea to check your output (including output of intermediate 
## steps) to make sure your code runs as you intended


### Summary table, by group ###

# create summary statistics
caffeine %>% 
  summarize(
    n = n(), 
    n_matcha = sum(matcha), 
    score_range = score %>% {max(.) - min(.)}
  )

# summary statistics subset by coffee treatment status
caffeine %>% 
  group_by(coffee) %>%
  summarize(
    n = n(), 
    n_matcha = sum(matcha), 
    score_range = score %>% {max(.) - min(.)}
  ) %>%
  ungroup() ## it is best practice to ungroup() at the end

# summary statistics subset by coffee x matcha treatment status
caffeine %>% 
  group_by(coffee, matcha) %>%
  summarize(
    n = n(), 
    mean_score = mean(score), 
    score_range = score %>% {max(.) - min(.)}
  ) %>%
  ungroup()


########################################
## Data visualization with `ggplot2`
########################################

mtcars %>% 
  ggplot(mapping = aes(x = hp, y = mpg)) + 
  geom_smooth(linetype = 6, se = FALSE) +
  geom_point(size = 3) + 
  facet_grid(
    am ~ cyl,
    labeller = labeller(
      am = c(`0` = "auto transmission", `1` = "manual transmission"), 
      cyl = c(`4` = "4 cylinders", `6` = "6 cylinders", `8` = "8 cylinders")
    )
  ) + 
  scale_x_continuous(breaks = c(100, 200, 300)) + 
  scale_y_continuous(limits = c(0, 40), n.breaks = 3) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        legend.position = "none") + 
  labs(x = "horsepower", y = "miles per gallon", 
       title = "Relationship between horsepower and fuel consumption", 
       subtitle = "by number of cylinders and transmission type")
