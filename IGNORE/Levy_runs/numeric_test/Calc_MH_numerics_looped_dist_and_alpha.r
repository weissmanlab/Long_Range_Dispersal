library(libstableR)
library(argparse)
library(mdatools)
args <- commandArgs(trailingOnly = TRUE)
if (length(args)!=4) {

  print(length(args))
  stop("Four input arguments are required -  number of timesteps, scale parameter, mutation rate, and rho inverse", call.=FALSE)
  
}
# store name of directory and change into said directory


#alpha <- as.double(args[1])
#initial_distance <- as.double(args[2])
#dir <- paste("./alpha_value_", alpha, "/distance_value_", initial_distance, sep = "")
#setwd(dir)
num_time_steps <- as.integer(args[1]) # number of times considered
scale_parameter <- as.double(args[2])
MU <- as.double(args[3])
rho_inverse <- as.double(args[4])
rho <- 1/rho_inverse

for(alpha_counter in 1:4)
{ alpha <- 1.05 + .2*alpha_counter
for(dist_counter in 0:13)
{ 
	initial_distance = exp(dist_counter)
	
	# round dist values so that they match simulation dists exactly
	if(dist_counter == 0){initial_distance = 1} 
	if(dist_counter == 1){initial_distance = 2.7} 
	if(dist_counter == 2){initial_distance = 7.39} 
	if(dist_counter == 3){initial_distance = 20.09} 
	if(dist_counter == 4){initial_distance = 54.60} 
	if(dist_counter == 5){initial_distance = 148.41} 
	if(dist_counter == 6){initial_distance = 403.43} 
	if(dist_counter == 7){initial_distance = 1096.63} 
	if(dist_counter == 8){initial_distance = 2980.96} 
	if(dist_counter == 9){initial_distance = 8103.08} 
	if(dist_counter == 10){initial_distance = 22026} 
	if(dist_counter == 11){initial_distance = 59874} 
	if(dist_counter == 12){initial_distance = 162755} 
	
	## add zero dist
	if(dist_counter == 13){initial_distance = 0} 


X <- initial_distance #more convenient name for distance variable
Laplace_Domain_Kernel_at_X <- 0
Laplace_Domain_Kernel_at_ZERO <- 0
timestep_size <- 1/10
#timestep_size <- 1
#We start with the continous time expression and account for discretized time and finite width coalescence zone.
for (time_step in 1:num_time_steps) { 
time <- time_step*timestep_size

pars <- c(alpha, 0, scale_parameter*(time)**(1/alpha), 0)

 Laplace_Domain_Kernel_at_X <- Laplace_Domain_Kernel_at_X +  stable_pdf(X , pars)*exp(-2*MU*time)*timestep_size
 Laplace_Domain_Kernel_at_ZERO <- Laplace_Domain_Kernel_at_ZERO +  stable_pdf(0 , pars)*exp(-2*MU*time)*timestep_size  

            # coalescence zone has width one in our simulations

    }


#Mean_Homozygosity <-  Laplace_Domain_Kernel_at_X/(rho)
#Mean_Homozygosity <- rho_inverse*( Laplace_Domain_Kernel_at_X -  (rho_inverse*Laplace_Domain_Kernel_at_X**2)/(1 + rho_inverse*Laplace_Domain_Kernel_at_ZERO))
Mean_Homozygosity <- ( Laplace_Domain_Kernel_at_X )/(rho + Laplace_Domain_Kernel_at_ZERO)
#output_file_name <- paste("NUMERIC_mean_homozygosity", alpha , "distance", initial_distance, "MU", MU, "rho_inverse_" , rho_inverse , ".txt", sep = "") 
output_file_name_plot <- paste("NUMERIC_mean_homozygosity", alpha , "distance", initial_distance, "MU", MU, "rho_inverse_" , rho_inverse , ".txt", sep = "") 
#sink(output_file_name)
#constrained_probability <- stable_cdf(initial_distance + .5, constrained_pars) -   stable_cdf(initial_distance - .5, constrained_pars)
#cat(paste(initial_distance, MU, Mean_Homozygosity, Mean_Homozygosity, Mean_Homozygosity, sep = " "))
#sink()

#setwd("../../MH_plots")
#setwd("/MH_plots")
sink(output_file_name_plot)
#constrained_probability <- stable_cdf(initial_distance + .5, constrained_pars) -   stable_cdf(initial_distance - .5, constrained_pars)
cat(paste(alpha, initial_distance, rho_inverse, MU ,Mean_Homozygosity, Mean_Homozygosity, Mean_Homozygosity, sep = " "))
sink()

}
}

