#####################
### ZOO800: HW 10 ###
#####################

# Author: Rebekkah LaBlue
# Concept: Linear Regression & Assumption Checking
# Due: November 10, 2025


### --- PACKAGES --- ###

library(tidyverse)
library(dplyr)
library(purrr)
library(magrittr)
library(car)
library(units)
library(stats)
library(ggplot2)
library(viridis)
library(ggfortify)
library(here)


### --- PROBLEM --- ###

### Ordinary least squares regression is a valuable tool, but one that comes with several important
# assumptions that should be checked. Violations of these assumptions can sometimes result in biased
# estimates of the regression parameters with consequences for the accuracy of predictions.


###################
### OBJECTIVE 1 ###
###################

### Find real data related to your research focus that fits the assumptions of ordinary least squares
# regression, estimate the regression parameters, and compare model predictions at two X values.


### OPTIONS: regress PA on land cover types and tmaxes?
# Scale predictor var/s but not response (not REQUIRED for only one covr, but makes interp easier,
# ie. coefficient response size for x for each 1 unit increase of y...)
# Spp data are det/nondet, ie.binomial and so not usable for this assignment
# use autoplot() for diagnostic plots, as well as histogram, to assess assumptions


### --- A --- ###
### Using either your own data or data that you find in an online database find two continuous variables
# that might reasonably be hypothesized to have a causal association (i.e., one variable is clearly the
# response, Y, and the other the predictor, X) and have sufficient numbers of paired observations (> 30).

covars_raw <- read.csv(here("wibba_covars_raw.csv"), header = TRUE)
colnames(covars_raw)
lm_df <- covars_raw %>%
  dplyr::select(pa_percent, lat, atlas_block) %>%
  filter(!is.na(pa_percent), !is.na(lat))

### All my data is at the "block"-level, ie. 25km^2 survey blocks (atlas_block) mapped over 
# the entire state of Wisconsin for the Wisconsin Breeding Bird Atlas (N = 3337)
# Covars to use: pa_percent (percent of total protected area per block), lat (latitude centroids of Atlas blocks)


### --- B --- ###
### Using lm(), fit a linear regression to these data.

lm1 <- lm(pa_percent ~ lat , data = lm_df)
summary(lm1)


### --- C --- ###
### Evaluate the model residuals for signs that regression assumptions are violated. You should evaluate
# at least three assumptions and for each one state to what extent you believe it is violated and how you
# know. You should plot figures and write your response as comments embedded in the code.

autoplot(lm1, which = 1:6)
# Residuals v Fitted: Assumption of linearity is violated (bend in line and clumping of points)
# Q-Q: Assumption of linearity violated (tails skewed, esp right tail)
# Scale-location: Assumption of homoschedasticity violated (bent line, point clumping indicate non-constant variance; 
# variance increases with fitted values)
# Cook's Distance: only 3 points are potentially influential observations 
# Residuals v Leverage, Cooks v Leverage: No highly influential points indicate no worrisome outliers

### Overall: Key assumptions of linearity, homoschedasticity, and normality violated
# Transformations likely necesssary to obtain better model fit


### --- D --- ###
### Generate predictions and associated prediction intervals for two X values: one at the median of X and
# the other at the 95th percentile of X. How do the prediction intervals differ?

x_median <- median(lm_df$lat)
x_95pct <- quantile(lm_df$lat, probs = 0.95)

lm_df2 <- data.frame(lat = c(x_median, x_95pct))
preds <- predict(lm1, newdata = lm_df2, interval = "prediction", level = 0.95)
bind_cols(lm_df2, as.data.frame(preds))

### Interp.:
# As latitude increases so does the percent protected area (pa) in each survey Atlas block (fit value > 0).
# My prediction intervals around the fitted value for both the median and 95th percentile values
# are very large/wide, hardly differ (upr - lwr = ~102.6, and contain negative values--all of which
# indicate the linear model is a poor fit for this data--echoing the violated regression assumptions. 
# The intervals being roughly the same width suggests values/uncertainty at tails/extremes 
# are similar and therefore stable, but it's clear that this model is a poor fit. Negative values
# in particular are problematic because you can't have 'negative' protected area as far as my 
# response variable is concerned.


###################
### OBJECTIVE 2 ###
###################

### I’ve said in class that the regression parameter estimates are fairly robust to modest deviations from
# normality. However, estimates of uncertainty are more sensitive. Evaluate whether this is true.

### --- A --- ###
### Generate linear regression data where the error in the response variable Y is not quite normally
# distributed (but still unimodal). A lognormal or negative binomial distribution should work. No
# error in X this time. 100 X, Y pairs should be good.

set.seed(333)

n <- 100 # n pairs
slope <-  -0.7 # set a true slope
intercept <- 3.3 # set a true intercept
X <- runif(n, min = -7, max = 7) # create x from a uniform distribution
epsilon <- rlnorm(n, meanlog = 0, sdlog = 0.6) # create skewed error from lognorm dist

Y <- intercept + slope * X + epsilon

ggplot(data.frame(X, Y), aes(x = X, y = Y)) + 
  geom_point() +
  geom_smooth(method = "lm")


### --- B --- ###
### Fit a linear regression to the data.

lm_b <- lm(Y ~ X)
summary(lm_b)
autoplot(lm_b, which = 1:6)

### --- C, D, E --- ###
### Repeat this process and keep track of the true and estimated slope and intercept.

Simulate_lm_b <- function(n = 100, slope = -0.7, intercept = 3.3) { # set base/true values for function to run
  X <- runif(n, min = -0.7, max = 0.7) # confg X
  epsilon <- rlnorm(n, meanlog = 0, sdlog = 0.6) # config error/epsilon 
  Y <- intercept + slope * X + epsilon # equation
  fit <- lm(Y ~ X) # fit model on predictor, response
  coef_fit <- coef(fit) # fit coeficients for retrieval (slope, int, X)
  preds <- predict(fit, newdata = data.frame(X, Y), interval = "prediction", level = 0.95) # create prediction interval
  interval <- mean(Y >= preds[, "lwr"] & Y <= preds[, "upr"]) # set bounds 
  tibble( # store 
    sims = sims,
    est_intercept = coef_fit[1],
    est_slope = coef_fit[2],
    interval = interval
  )
}

set.seed(333)
lm_b_sims <- replicate(100, Simulate_lm_b(), simplify = FALSE) %>% # replicate function over 100 sims
  bind_rows()

sims_summary <- lm_b_sims %>% # extract means, sds for coefficients of interest 
  summarise(
    mean_intercept = mean(est_intercept),
    sd_intercept = sd(est_intercept),
    mean_slope = mean(est_slope),
    sd_slope = sd(est_slope),
    mean_interval = mean(interval)
  )

sims_summary # compare
# Pred slop = -0.684 +/- 0.195, True slope = -0.7 
# Pred intercept = 4.5 +/- 0.09, True intercept = 3.3
# Match: Pretty close!


### --- F --- ###
### What fraction of your data (Y values) falls within the 95% prediction interval?
# Roughly 95%! (0.954)


### --- G --- ###
### What does this imply for how your estimated uncertainty in the predictions compares to the
# true uncertainty.

# Overall, they are pretty close to one another with more simulations (CLT)
