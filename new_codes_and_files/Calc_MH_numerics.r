library(libstableR)
library(argparse)
library(mdatools)
args <- commandArgs(trailingOnly = TRUE)
if (length(args)!=6) {

  print(length(args))
  stop("Six input arguments are required - alpha, initial distance, number of timesteps, scale parameter, mutation rate, and rho inverse", call.=FALSE)
  
}
# store name of directory and change into said directory


alpha <- as.double(args[1])
initial_distance <- as.double(args[2])
dir <- paste("./alpha_value_", alpha, "/distance_value_", initial_distance, sep = "")
setwd(dir)
num_time_steps <- as.integer(args[3]) # number of times considered
scale_parameter <- as.double(args[4])
MU <- as.double(args[5])
rho_inverse <- as.double(args[6])
X <- initial_distance #more convenient name for distance variable
Laplace_Domain_Kernel_at_X <- 0
Laplace_Domain_Kernel_at_ZERO <- 0
timestep_size <- 1/1000
#We start with the continous time expression and account for discretized time and finite width coalescence zone.
for (time in 1:(1/timestep_size)*(num_time_steps)) { 
pars <- c(alpha, 0, scale_parameter*(time*timestep_size)**(1/alpha), 0)

 Laplace_Domain_Kernel_at_X <- Laplace_Domain_Kernel_at_X +  stable_pdf(X , pars)*exp(-2*MU*time*timestep_size)*timestep_size
 Laplace_Domain_Kernel_at_ZERO <- Laplace_Domain_Kernel_at_ZERO +  stable_pdf(0 , pars)*exp(-2*MU*time*timestep_size)*timestep_size  

            # coalescence zone has width one in our simulations

    }


Mean_Homozygosity <- rho_inverse*( Laplace_Domain_Kernel_at_X -  (rho_inverse*Laplace_Domain_Kernel_at_X**2)/(1 + rho_inverse*Laplace_Domain_Kernel_at_ZERO))
output_file_name <- paste("NUMERIC_mean_homozygosity", alpha , "distance", initial_distance, "MU", MU, "rho_inverse_" , rho_inverse , ".txt", sep = "") 
output_file_name_plot <- paste("NUMERIC_mean_homozygosity", alpha , "distance", initial_distance, "MU", MU, "rho_inverse_" , rho_inverse , ".txt", sep = "") 
sink(output_file_name)
#constrained_probability <- stable_cdf(initial_distance + .5, constrained_pars) -   stable_cdf(initial_distance - .5, constrained_pars)
cat(paste(initial_distance, MU, Mean_Homozygosity, Mean_Homozygosity, Mean_Homozygosity, sep = " "))
sink()

setwd("../../MH_plots")

sink(output_file_name_plot)
#constrained_probability <- stable_cdf(initial_distance + .5, constrained_pars) -   stable_cdf(initial_distance - .5, constrained_pars)
cat(paste(alpha,   MU , rho_inverse,   initial_distance, Mean_Homozygosity, Mean_Homozygosity, Mean_Homozygosity, sep = " "))
sink()



setwd("..")
