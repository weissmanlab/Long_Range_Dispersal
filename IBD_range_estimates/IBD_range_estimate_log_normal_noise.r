library(tidyverse)
library(modelr)
options(na.action = na.warn)
power_law_model  <- function(x, xbar, rhomu, alpha) {(gamma(abs(alpha) + 1)/(2*pi))*sin(pi*abs(alpha)/2)*((x/abs(xbar))^(-abs(alpha) -1))/abs(xbar*rhomu)}
exponential_model  <- function(x, xbar, rhomu){exp(-x/abs(xbar))/(4*abs(rhomu*xbar) + 1)}
num_ending_distances <- 4
num_sigma <- 5
total_length <- num_sigma*num_ending_distances
starting_distances <- rep(1, total_length)
ending_distances <- c(2*rep(1,num_sigma), 10*rep(1,num_sigma), 20*rep(1,num_sigma), 40*rep(1,num_sigma))
sigma_values <- rep(.02*2^c(1:num_sigma), num_ending_distances)
zerovec <- 0*c(1:total_length)
num_dist_points <- 5

LLR_df <- tibble(starting_distance = starting_distances, ending_distance = ending_distances, sigma =  sigma_values, LLR_avg =  zerovec, LLR_low =  zerovec, LLR_high =  zerovec)

for(i in 1:total_length){
  start_dist <- LLR_df$starting_distance[i]
  end_dist <- LLR_df$ending_distance[i]
  distance_interval <- end_dist - start_dist
  sigma_val <- LLR_df$sigma[i]
  dist_vec <- start_dist*rep(1,num_dist_points) + (distance_interval/(num_dist_points-1))*c(0:(num_dist_points-1))
  means <- power_law_model(dist_vec, 1, 1, 1.5)
  MLE_power_avg <- 0
  MLE_exp_avg <- 0
  num_rand_trials <- 500
  for(Q in 1:num_rand_trials){
  simulated_data <- rlnorm(num_dist_points, meanlog = log(means), sdlog = sigma_val)
  
  neg_power_LL_func <- function(par) # par is c(xbar, rhomu, alpha, SIGMA) # bound alpha at 2
   {-1*sum(dlnorm(simulated_data, meanlog =log(power_law_model(dist_vec, par[1], par[2], 2/(1 + exp(par[3]))  )), sdlog = abs(par[4]), log = TRUE))}
  
  neg_exponential_LL_func <- function(par) #par is c(xbar, rhomu, SIGMA)
   {-1*sum(dlnorm(simulated_data, meanlog =log(exponential_model(dist_vec, par[1], par[2])), sdlog = abs(par[3]), log = TRUE))}
  
  MLE_power <- -1*optim(c(1, 1, 1.5, 1), neg_power_LL_func)$value
  MLE_exp <- -1*optim(c(5, 1, 1), neg_exponential_LL_func)$value
  MLE_power_avg <- MLE_power_avg + MLE_power/num_rand_trials
  MLE_exp_avg <- MLE_exp_avg + MLE_exp/num_rand_trials
  }
  
  Log_of_Likelihood_ratio <- MLE_power_avg - MLE_exp_avg
  LLR_df$LLR_avg[i] <- Log_of_Likelihood_ratio
  }
LLR_df <- mutate(LLR_df, sigma = as.factor(sigma))
ggplot(data = LLR_df) + geom_line(aes(x = ending_distance, y = LLR_avg, color = sigma)) + labs(x= "Size of Spatial Range", y = "Log of Likelihood Ratio for IBD Models" )

