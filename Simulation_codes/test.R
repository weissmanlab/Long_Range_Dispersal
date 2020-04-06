library(likelihood)
library(ggplot2)
library(stats)
require(lattice)
require(reshape2)
library(pracma)
library(Deriv)
library(tidyverse)
test_data <- as_tibble((read.table("aggregated_mh_data.txt",  header = T)))
test_data  <- mutate(test_data, xbar = c*(2*mu)^(-1/alpha), x_nat = x0/xbar )
test_data  <- mutate(test_data, rho = (coal)^-1 )


test_data  <- filter(test_data, alpha == 1.22,  mu < 1, c ==20, x_nat < 40)
test_data  <- mutate(test_data, hom_scale_approx = (rho*xbar*mu ), Shom =hom_scale_approx*hom, ShomLow =hom_scale_approx*homLow, ShomHigh =hom_scale_approx*homHigh  )
#



Levy_model_1D_asymptotic	<- function(x, xbar, rhomu, alpha) {
  delta <- 10^-3*abs(xbar)
  
  if(alpha > 1) {psi0 =  (2*alpha*sin(pi/alpha)*abs(rhomu*xbar) +1)^-1}
  if(alpha == 1){psi0 = log(abs(xbar)/delta)/(log(abs(xbar)/delta) + 2*pi*abs(rhomu*xbar)) }
  if(alpha < 1) {psi0 = 2^((1-alpha)/2)*gamma(.5 - alpha/2)/( 2^((1-alpha)/2)*gamma(.5 - alpha/2) + 4*pi*abs(rhomu)*abs(xbar)^alpha*delta^(1-alpha))  }
  asymptotic_series <- 0 
  for(Q in 1:1){	 	
    asymptotic_series <- asymptotic_series + (-1)^(Q -1)*(2*pi*rhomu*xbar)^-1*(1 - psi0)*sin(Q*alpha*pi/2)*gamma(1 + Q*alpha)*abs(x/xbar)^(-Q*alpha -1)}
  return(asymptotic_series)
}



# get psi0 at mu = 0.  This is


test_data  <- mutate(test_data, Levy_model =Levy_model_1D_asymptotic(x0, xbar,rho*mu, alpha), SLevy_model = hom_scale_approx*Levy_model    )

test_data  <- mutate(test_data,mu = as.factor(mu) )



ggplot(data = test_data) + geom_pointrange(aes(x = x_nat, y = Shom, ymin=ShomLow, ymax = ShomHigh, color = mu)) + 
  scale_y_log10() + scale_x_log10() + geom_point(aes(x = x_nat, y = SLevy_model))

test_data_CDF <- as_tibble((read.table("aggregated_cdf_data.txt",  header = T)))
test_data_CDF <- filter(test_data_CDF, alpha < 1)

# ALPHA LESS THAN ONE ONLY. FOR ALPHA >= 1 max_CDF = 1
max_CDF <- function( alpha, rho, Da) {
 
  delta <- 1
  psi0 = 2^((1-alpha)/2)*gamma(.5 - alpha/2)/( 2^((1-alpha)/2)*gamma(.5 - alpha/2) + 4*pi*abs(rho)*abs(Da)*delta^(1-alpha))  
  return(psi0)
  
  
  }


test_data_CDF <- mutate(test_data_CDF, rho = coal^-1, Da = .5*c^alpha, 
                        comp_cdf = max_CDF( alpha, rho, Da) - cdf,
                       comp_cdfLow = max_CDF( alpha, rho, Da) - cdfHigh,
                       comp_cdfHigh = max_CDF( alpha, rho, Da) - cdfLow, 
                        Scaled_T  = (2*(rho^alpha)*Da)^((1-alpha)^-1)*T, 
                        Scaled_comp_cdf = pi*(1/alpha - 1)/(gamma(1/alpha + 1))*comp_cdf,
                        Scaled_comp_cdfLow = pi*(1/alpha - 1)/(gamma(1/alpha + 1))*comp_cdfLow,
                        Scaled_comp_cdfHigh = pi*(1/alpha - 1)/(gamma(1/alpha + 1))*comp_cdfHigh,
                        model = Scaled_T^(1-1/alpha)  )

#CDF is improper for alpha less than one!  We can't just do 1 - cdf

test_data_CDF <- mutate(test_data_CDF, alpha =as.factor(alpha))
ggplot(data = test_data_CDF) + geom_pointrange(aes(x = Scaled_T, y = Scaled_comp_cdf, ymin = Scaled_comp_cdfLow, ymax = Scaled_comp_cdfHigh, color = alpha )) + 
scale_y_log10() + scale_x_log10() +geom_line(aes(x = Scaled_T, y =model, color = alpha))


# Good agreement for such a small number of trials


