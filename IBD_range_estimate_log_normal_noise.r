library(ggplot2)
require(reshape2)

# To generate noisy data, we use a 1D power law model with alpha = 1.5, xbar = 1, mu = .1, and rho = 10
powerfunction <- function(x){(gamma(2.5)/(2*pi))*sin(pi*1.5/2)*x^-2.5}
expfunction <- function(x){exp(-x)/(4*10*.1*1 + 1)}
num_sigma <- 8 # number of distinct noise strengths considered.  We want to see how differences in model log likelihoods vary with the amount of noise
index_to_sigma <- data.frame(1, num_sigma)   # relates integer index to noise in subset of data
num_distance_points <- 5 # number of distinct spatial points in sample
num_point_spreads <- 4   # points are spread evenly across log distances.  
# each increment increases the point spread by half an order of magnitude, 10^.5
starting_point <- 1
# all sets of points start at starting_point
for(i in 1:num_sigma){       # defines index to sigma
if(i == 1){index_to_sigma[1, i] = 0} # dummy entry - this column has distance axis
if(i == 2){index_to_sigma[1, i] = .01}
if(i == 3){index_to_sigma[1, i] = .05}
if(i == 4){index_to_sigma[1, i] = .1}
if(i == 5){index_to_sigma[1, i] = .25}
if(i == 6){index_to_sigma[1, i] = .5}
if(i == 7){index_to_sigma[1, i] = 1}
if(i == 8){index_to_sigma[1, i] = 2.5}
if(i > 8){index_to_sigma[1, i] = 5}
}
#power_law_model_linear_input  <- function(x, xbar, rho, mu) {(gamma(2.5)/(2*pi))*sin(pi*1.5/2)*((x/abs(xbar))^-2.5)/abs(xbar*rho*mu)}
power_law_model_linear_input  <- function(x, xbar, rho, mu, alpha) {(gamma(abs(alpha) + 1)/(2*pi))*sin(pi*abs(alpha)/2)*((x/abs(xbar))^(-abs(alpha) -1))/abs(xbar*rho*mu)}
exponential_model_linear_input  <- function(x, xbar, rho, mu){exp(-x/abs(xbar))/(4*abs(rho*mu*xbar) + 1)}
#power_law_model_log_input  <- function(logx, xbar, rho, mu) {(gamma(2.5)/(2*pi))*sin(pi*1.5/2)*((exp(logx)/abs(xbar))^-2.5)/abs(xbar*rho*mu)}
power_law_model_log_input  <- function(logx, xbar, rho, mu, alpha) {(gamma(abs(alpha) + 1)/(2*pi))*sin(pi*abs(alpha)/2)*((exp(logx)/abs(xbar))^(-abs(alpha) -1))/abs(xbar*rho*mu)}
exponential_model_log_input  <- function(logx, xbar, rho, mu){exp(-exp(logx)/abs(xbar))/(4*abs(rho*mu*xbar) + 1)}
relative_log_likelihood_log_spread_avg <- data.frame(num_point_spreads, num_sigma)
for(i in 1:num_sigma){ for(j in 1:num_point_spreads){relative_log_likelihood_log_spread_avg[j, i] <- 0 }}
# we want to see how the difference in log likelihood between the models varies with the noise in the data and the spread between a fixed number of points.
Num_trials = 100 # We can average over many different random samples to get a sense of expected discrepancy betweeen model fits






for(Q in 1: Num_trials){
    relative_log_likelihood_log_spread <- data.frame(num_point_spreads, num_sigma)
     for( g in 1:num_point_spreads){ 
       relative_log_likelihood_log_spread[g, 1] <- starting_point*exp((g)*.5*log(10))  # this column contains the end point
       power_func_w_noise_log_spread <- data.frame( num_distance_points, num_sigma)
       for(i in 2:num_sigma){ for(j in 1:num_distance_points){
       power_func_w_noise_log_spread[j, 1] <- j # this column just contains log distance
        sigma_var <- index_to_sigma[1,i]
          power_func_w_noise_log_spread[j, i] <- rlnorm(1, meanlog = log(powerfunction(starting_point*exp((g/num_distance_points)*j*.5*log(10)) )) - sigma_var^2/2, sdlog = sigma_var) # this randomly generates noisy data at multiple distances
          }}
    for(i in 2:num_sigma){ 
     dummy_func_power <-  function( xbar, rho, mu, alpha, sigma) {sum((log(power_law_model_log_input( log(starting_point)+(g/num_distance_points)*.5*log(10)*(1:num_distance_points), xbar, rho, mu, alpha)) - (sigma^2/2) -log(power_func_w_noise_log_spread[ 1:num_distance_points, i]))^2/(2*( sigma^2)) + log((2*pi)^.5*abs(sigma)*power_func_w_noise_log_spread[ 1:num_distance_points, i])   ) }
     dummy_func_exponential <-  function( xbar, rho, mu, sigma) {sum((log(exponential_model_log_input( log(starting_point)+(g/num_distance_points)*.5*log(10)*(1:num_distance_points), xbar, rho, mu)) -log(power_func_w_noise_log_spread[ 1:num_distance_points, i]))^2/(2*( sigma^2)) + log((2*pi)^.5*abs(sigma)*power_func_w_noise_log_spread[ 1:num_distance_points, i])   ) }
     dummy_func_power_par <- function(par){ dummy_func_power(par[1], par[2], par[3], par[4],  par[5]) }
     dummy_func_power_par_unconstrained <- function(par) { dummy_func_power(par[1], par[2], par[3], 2/(1 + exp(par[4])), par[5]) }
     # we've reparameterized, replacing alpha with 2*logistic(x).  Now we can do unconstrained optimization
     dummy_func_exponential_par <- function(par){ dummy_func_exponential(par[1], par[2], par[3], par[4]) }
     
     power_log_likelihood_log_spread <- -1*optim(par =c(10, 1, 1, 1, 1), dummy_func_power_par_unconstrained )$value  
     power_AIC_log_spread <-  2*5 -  2*power_log_likelihood_log_spread  
     exponential_log_likelihood_log_spread <- -1*optim(par =c(10, 1, 1, 1), dummy_func_exponential_par )$value
     exponential_AIC_log_spread <-  2*4 -  2*exponential_log_likelihood_log_spread  
     
     relative_log_likelihood_log_spread[g, i] <-  (exponential_AIC_log_spread - power_AIC_log_spread)/2 
     
      }}

    relative_log_likelihood_log_spread_avg <- relative_log_likelihood_log_spread_avg +  relative_log_likelihood_log_spread/Num_trials
     

}



names(relative_log_likelihood_log_spread_avg)[1] <- "endpoint"
names(relative_log_likelihood_log_spread_avg)[2] <- "sigma .01"
names(relative_log_likelihood_log_spread_avg)[3] <- "sigma .05"
names(relative_log_likelihood_log_spread_avg)[4] <- "sigma .1"
names(relative_log_likelihood_log_spread_avg)[5] <- "sigma .25"
names(relative_log_likelihood_log_spread_avg)[6] <- "sigma .5"
names(relative_log_likelihood_log_spread_avg)[7] <- "sigma 1"
names(relative_log_likelihood_log_spread_avg)[8] <- "sigma 2.5"
#df <- data.frame(time = 1:10,
#                 a = cumsum(rnorm(10)),
#                 b = cumsum(rnorm(10)),
#                 c = cumsum(rnorm(10)))
#df <- melt(df ,  id.vars = 'time', variable.name = 'series')

df <- melt(relative_log_likelihood_log_spread_avg,  id.vars = 'endpoint', variable.name = 'series')
print(df)
ggplot(df, aes(endpoint,value)) + geom_line(aes(colour = series))
#print(colnames(relative_log_likelihood_log_spread_avg))
#extend array to more values of noise and distance, add linear spread array in addition to log spread array
#plot(as.vector(relative_log_likelihood_log_spread_avg[, 6]))
#add error bars for noise in the log likelihod for a given sample?


   