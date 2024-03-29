library(tidyverse)
library(modelr)
options(na.action = na.warn)

#First we define the two IBD models that we'll be comparing

power_law_model  <- function(x, xbar, rhomu, alpha) {(gamma(abs(alpha) + 1)/(2*pi))*sin(pi*abs(alpha)/2)*((x/abs(xbar))^(-abs(alpha) -1))/abs(xbar*rhomu)}
exponential_model  <- function(x, xbar, rhomu){exp(-x/abs(xbar))/(4*abs(rhomu*xbar) + 1)}


# We will compare the fits of the two models to simulated data drawn from the power law model + noise
# We want see how the likeilhood ratio of the two models varies with the amound of noise and the size of the spatial range considered

num_ending_distances <- 4   # number of different range sizes
num_sigma <- 4 #number of different noise sizes
total_length <- num_sigma*num_ending_distances #total number of rows of dataframe
starting_distances <- rep(1, total_length) # create vector with one starting distance, x = 1
ending_distances <- c(2*rep(1,num_sigma), 10*rep(1,num_sigma), 20*rep(1,num_sigma), 40*rep(1,num_sigma)) # create vector with all ending distances
sigma_values <- rep(.02*2^c(1:num_sigma), num_ending_distances) # create vector with all noise values
zerovec <- 0*c(1:total_length)
num_dist_points <- 5   # number of actual points in the distance interval

# Next we aggregate all of these parameter values into a dataframe and calculate the log of likelihood ratio for each choice of parameters
LLR_df <- tibble(starting_distance = starting_distances, ending_distance = ending_distances, sigma =  sigma_values, LLR_avg =  zerovec, LLR_low =  zerovec, LLR_high =  zerovec)


for(i in 1:total_length){  # find likelihood of models for each row in LLR_df
  start_dist <- LLR_df$starting_distance[i]
  end_dist <- LLR_df$ending_distance[i]
  distance_interval <- end_dist - start_dist
  sigma_val <- LLR_df$sigma[i]
  dist_vec <- start_dist*rep(1,num_dist_points) + (distance_interval/(num_dist_points-1))*c(0:(num_dist_points-1))
  means <- power_law_model(dist_vec, 1, 1, 1.5)
  MLE_power_avg <- 0
  MLE_exp_avg <- 0
  num_rand_trials <- 500
  Trials_df <- tibble(MLE_power = 0*c(1:num_rand_trials), MLE_exp= 0*c(1:num_rand_trials), LLR= 0*c(1:num_rand_trials))
  for(Q in 1:num_rand_trials){  # Average over many trials of noisy data
  simulated_data <- rlnorm(num_dist_points, meanlog = log(means), sdlog = sigma_val) # We generate random data, power law model + noise
  
  neg_power_LL_func <- function(par) #  Log likeihood for power law model on noisy data
   {-1*sum(dlnorm(simulated_data, meanlog =log(power_law_model(dist_vec, par[1], par[2], 2/(1 + exp(par[3]))  )), sdlog = abs(par[4]), log = TRUE))}
  
  neg_exponential_LL_func <- function(par) #  Log likeihood for exponential model on noisy data
   {-1*sum(dlnorm(simulated_data, meanlog =log(exponential_model(dist_vec, par[1], par[2])), sdlog = abs(par[3]), log = TRUE))}
  
  Trials_df$MLE_power[Q] <- -1*optim(c(1, 1, 1.5, 1), neg_power_LL_func)$value # maximum likelihood estimate for power law model on data from this trial
  Trials_df$MLE_exp[Q] <- -1*optim(c(5, 1, 1), neg_exponential_LL_func)$value # maximum likelihood estimate for exponential model on data from this trial
  Trials_df$LLR[Q] <- Trials_df$MLE_power[Q] - Trials_df$MLE_exp[Q] 
  }
  
  LLR_df$LLR_avg[i] <- mean(Trials_df$LLR)  # We have log of likelihood ratio averaged over many noisy data sets for each choice of noise size and range size
  LLR_df$LLR_low[i] <- mean(Trials_df$LLR) - sd(Trials_df$LLR)
  LLR_df$LLR_high[i] <- mean(Trials_df$LLR) + sd(Trials_df$LLR)
  }
LLR_df <- mutate(LLR_df, sigma = as.factor(sigma))
ggplot(data = LLR_df) + geom_pointrange(aes(x = ending_distance, y = LLR_avg, ymin = LLR_low, ymax = LLR_high, color = sigma))  +
  geom_line(aes(x = ending_distance, y = LLR_avg,  color = sigma)) +
  labs(x= "Size of Spatial Range", y = "Log of Likelihood Ratio for IBD Models" ) +facet_wrap(~sigma)

