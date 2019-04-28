#########################################################
#install.packages("argparse")
#install.packages("libstableR")
#install.packages("mdatools")
library(libstableR)
library(argparse)
library(mdatools)
Sys.umask('000')



args <- commandArgs(trailingOnly = TRUE)
if (length(args)!=6) {

  print(length(args))
  stop("Six input arguments are required - alpha, initial distance, number of trials , number of timesteps, scale parameter, and mutation rate", call.=FALSE)
  
}

# store name of directory and change into said directory



alpha <- as.double(args[1])
initial_distance <- as.double(args[2])
dir <- paste("./alpha_value_", alpha, "/distance_value_", initial_distance, sep = "")
setwd(dir)
#system(paste("cd", dir, sep = " "))
num_trials <- as.integer(args[3]) # number of trials for each time
num_time_steps <- as.integer(args[4]) # number of times considered
weight_list <- data.frame(matrix(0, ncol = 1, nrow = num_trials))
scale_parameter <- as.double(args[5])
MU <- as.double(args[6])
pars <- c(alpha, 0, scale_parameter, 0) # These are the parameters for the step distribution of our discretized Levy flight

## weights remove bias from individual trajectories

for (trial in 1:num_trials) {
   input_file_name <- paste("Extra_step_jump_sizealpha", alpha , "distance", initial_distance, "MU", MU, "trial", (trial -1), ".txt", sep = "") 
   file <- read.table(input_file_name)
   
   extra_step <- file[1,1]
   #extra_step <- file[2, 1]
   #free_step <- file[1,1]
   #constrained_step <- file[2,1]
   #probability_of_origin_time_and_big_jump_time <-  file[3,1]
   #fake_probability <- stable_pdf(free_step, pars)*probability_of_origin_time_and_big_jump_time  # the dist used to sample the important paths 
   #true_probability <- stable_pdf(constrained_step, pars) # the prob under dist we really want to draw from
   ## ALL the other steps in the trajectories are drawn with the same probability in our true and fake distributions, so we only have to account for the bias in this one step
   weight <- stable_pdf(extra_step, pars)/1#dunif(extra_step, extra_step_mean -.5, extra_step_mean + .5) #(stable_cdf(extra_step+.5, pars) - stable_cdf(extra_step-.5, pars)) /1 # true_probability of entering coalescence zone /fake_probability  # this ratio corrects for the bias in our importance sampling procedure.
    weight_list[trial, 1] <- weight
  file.remove(input_file_name) #free up disk space
   }
 
#weight_sum <- sum(weight_list[,1])  DO NOT NORMALIZE.  UN-NORMALIZED weights w_i/num_trials are what give us the proper estimate of the expectation over all paths
#weight_list[,1] <- weight_list[,1]/weight_sum #normalize list




for (trial in 1:num_trials) {
    
   weight <- weight_list[trial, 1]
   output_file_name <- paste("trajectory_WEIGHT_alpha", alpha , "distance", initial_distance, "MU", MU, "trial", (trial -1), ".txt", sep = "") 
   sink(output_file_name)
   cat(weight)
   sink()

    
   }


## Above weights remove bias from individual trajectories

## Below calculates probability of path hitting coalescence zone.  
## Use this to convert from conditional distribution of origin-hitting paths to marginal distribution of all paths



##no_hit_probability <- 1 # probability of an arbitrary never hitting the coalescence zone num_time_steps time.  
##timestep_size <- 1/1000 ## dt
#We use the continous time approximation.
##for (time in 1:(1/timestep_size)*(num_time_steps )){ 
##constrained_pars <- c(alpha, 0, scale_parameter*(time*timestep_size)**(1/alpha), 0)

 # no_hit_probability <-  no_hit_probability *(1 - stable_cdf(initial_distance + .5, constrained_pars) -   stable_cdf(initial_distance - .5, constrained_pars)*timestep_size )
            # coalescence zone has width one in our simulations

 #   }


## NO EXTRA PROBABILITY FACTOR NEEDED. UN_NORMALIZED WEIGHTS ALREADY SAMPLE THE FULL DISTRIBUTION OF ALL PATHS
#probability_of_hitting_coalescence_zone <- 1 - no_hit_probability
#output_file_name <- paste("probability_of_constrained_trajectory_WEIGHT_alpha", alpha , "distance", initial_distance,  ".txt", sep = "") 
#sink(output_file_name)
#cat(probability_of_hitting_coalescence_zone)
#sink()



#system("cd ../..")
setwd("../..")
