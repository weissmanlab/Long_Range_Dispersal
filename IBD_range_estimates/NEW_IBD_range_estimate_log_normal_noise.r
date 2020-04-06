library(tidyverse)
require(reshape2)

# To generate noisy data, we use a 1D power law model 
#with alpha = 1.5, xbar = 1, mu = .1, and rho = 10

power_law_model  <- function(x, xbar, rho, mu, alpha) {(gamma(abs(alpha) + 1)/(2*pi))*sin(pi*abs(alpha)/2)*((x/abs(xbar))^(-abs(alpha) -1))/abs(xbar*rho*mu)}
exponential_model  <- function(x, xbar, rho, mu){exp(-x/abs(xbar))/(4*abs(rho*mu*xbar) + 1)}

