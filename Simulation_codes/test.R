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

test_data  <- mutate(test_data, Levy_model =Levy_model_1D_asymptotic(x0, xbar,rho*mu, alpha), SLevy_model = hom_scale_approx*Levy_model    )

test_data  <- mutate(test_data,mu = as.factor(mu) )



ggplot(data = test_data) + geom_pointrange(aes(x = x_nat, y = Shom, ymin=ShomLow, ymax = ShomHigh, color = mu)) + 
  scale_y_log10() + scale_x_log10() + geom_point(aes(x = x_nat, y = SLevy_model))

