library(tidyverse)
test_data <- as_tibble((read.table("small_dispersal_data.txt",  header = T)))
test_data  <- mutate(test_data, xbar = c*(2*mu)^(-1/alpha), 
                     rho = (coal)^-1, 
                     Da = c^alpha/2, 
                     ybar1 = Da^(1/alpha) )

Levy_model_1D_finite_dist	<- function(x, Da, rho, alpha) {
    model <-  gamma(1-alpha)*sin(pi*alpha/2)*x^(alpha - 1)/(2*pi*rho*Da)
}

Levy_model_1D_zero_dist	<- function(Da, rho, alpha) {
  delta <- 1
  psi0 <- 2^((1-alpha)/2)*gamma(.5 - alpha/2)/( 2^((1-alpha)/2)*gamma(.5 - alpha/2) + 4*pi*abs(rho)*abs(Da)*delta^(1-alpha))  
  return(psi0)
  }

#add columns to tibble
test_data  <- mutate(test_data, finite_dist_model =Levy_model_1D_finite_dist(x0, Da, rho, alpha), zero_dist_model = Levy_model_1D_zero_dist(Da, rho, alpha),mu = as.factor(mu))

test_data <- filter(test_data, mu != 1)

ggplot(data = test_data) + geom_pointrange(aes(x = x0, y = hom, ymin= homLow, ymax = homHigh, color = mu)) + 
  scale_y_log10() + scale_x_log10() + geom_line(aes(x = x0, y = finite_dist_model)) + geom_line(aes(x = x0, y = zero_dist_model), color = "blue")

