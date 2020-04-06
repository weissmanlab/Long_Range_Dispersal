library(tidyverse)
test_data <- as_tibble((read.table("aggregated_mh_data.txt",  header = T)))
test_data  <- mutate(test_data, xbar = c*(2*mu)^(-1/alpha), x_nat = x0/xbar, rho = (coal)^-1 )
test_data  <- filter(test_data, alpha == 1.22,  mu < 1, c ==20, x_nat < 40)
test_data  <- mutate(test_data, hom_scale_approx = (rho*xbar*mu ), Scaled_hom =hom_scale_approx*hom, Scaled_homLow =hom_scale_approx*homLow, Scaled_homHigh =hom_scale_approx*homHigh  )

Levy_model_1D_asymptotic	<- function(x, xbar, rhomu, alpha) {
    model <-  (2*pi*rhomu*xbar)^-1*sin(alpha*pi/2)*gamma(1 + alpha)*abs(x/xbar)^(-alpha -1)
}

#add columns to tibble
test_data  <- mutate(test_data, Levy_model =Levy_model_1D_asymptotic(x0, xbar,rho*mu, alpha), SLevy_model = hom_scale_approx*Levy_model,mu = as.factor(mu))

ggplot(data = test_data) + geom_pointrange(aes(x = x_nat, y = Scaled_hom, ymin=Scaled_homLow, ymax = Scaled_homHigh, color = mu)) + 
  scale_y_log10() + scale_x_log10() + geom_line(aes(x = x_nat, y = SLevy_model))

### ABOVE IS MH ANALYSIS.  BELOW IS CDF ANALYSIS.######################################

test_data_CDF <- as_tibble((read.table("aggregated_cdf_data.txt",  header = T)))
test_data_CDF <- filter(test_data_CDF, alpha < 1)

max_CDF <- function( alpha, rho, Da) { # ALPHA LESS THAN ONE ONLY. For alpha >= 1, max_CDF = 1.  Distribution of coalescence times is improper for alpha < 1.
  delta <- 1
  psi0 = 2^((1-alpha)/2)*gamma(.5 - alpha/2)/( 2^((1-alpha)/2)*gamma(.5 - alpha/2) + 4*pi*abs(rho)*abs(Da)*delta^(1-alpha))  
  return(psi0)
  }
## add new columns to the tibble
test_data_CDF <- mutate(test_data_CDF, rho = coal^-1, Da = .5*c^alpha, 
                        comp_cdf = max_CDF( alpha, rho, Da) - cdf,
                        comp_cdfLow = max_CDF( alpha, rho, Da) - cdfHigh,
                        comp_cdfHigh = max_CDF( alpha, rho, Da) - cdfLow, 
                        Scaled_T  = (2*(rho^alpha)*Da)^((1-alpha)^-1)*T, 
                        Scaled_comp_cdf = pi*(1/alpha - 1)/(gamma(1/alpha + 1))*comp_cdf,
                        Scaled_comp_cdfLow = pi*(1/alpha - 1)/(gamma(1/alpha + 1))*comp_cdfLow,
                        Scaled_comp_cdfHigh = pi*(1/alpha - 1)/(gamma(1/alpha + 1))*comp_cdfHigh,
                        model = Scaled_T^(1-1/alpha) , alpha =as.factor(alpha) )

#ggplot(data = test_data_CDF) + geom_pointrange(aes(x = Scaled_T, y = Scaled_comp_cdf, ymin = Scaled_comp_cdfLow, ymax = Scaled_comp_cdfHigh, color = alpha )) + 
#scale_y_log10() + scale_x_log10() +geom_line(aes(x = Scaled_T, y =model, color = alpha))



