#########################################################
#install.packages("argparse")
library(libstableR)
library(argparse)
library(mdatools)
args <- commandArgs(trailingOnly = TRUE)
if (length(args)!=6) {

  print(length(args))
  stop("Six input arguments are required - alpha, initial distance, number of trials , number of timesteps, scale parameter, and mutation rate", call.=FALSE)
  
}


alpha <- as.double(args[1])
initial_distance <- as.double(args[2])
num_trials <- as.integer(args[3]) # number of trials for each time
num_time_steps <- as.integer(args[4]) # number of times considered
weight_list <- data.frame(matrix(0, ncol = 1, nrow = num_trials))






scale_parameter <- as.double(args[5])
MU <- as.double(args[6])
pars <- c(alpha, 0, scale_parameter, 0)
# store name of directory and change into said directory
dir <- paste("./alpha_value_", alpha, "/distance_value_", initial_distance, sep = "")
#print(dir)
setwd(dir)
#system(cd_to_dir)

for (trial in 1:num_trials) {
    
   filestring <- paste("free_and_contrained_CHOSEN_TIME_jump_sizesalpha", alpha , "distance", initial_distance, "MU", MU, "trial", (trial -1), ".txt", sep = "") 
   #print(filestring)
   # -1 adjusts for fact that c++ loops start at zero and r loops start at 1
   file <- read.table(filestring)
   free_step <- file[1,1]
   constrained_step <- file[2,1]
   
   weight <- stable_pdf(constrained_step, pars)/stable_pdf(free_step , pars)
   weight_list[trial, 1] <- weight
   
   
    
   }


weight_sum <- sum(weight_list[,1])
weight_list[,1] <- weight_list[,1]/weight_sum #normalize lise



no_hit_probability <- 1 # probability of an arbitrary never hitting the coalescence zone num_time_steps time.  
timestep_size <- 1/1000
#We start with the continous time expression and account for discretized time and finite width coalescence zone.
for (time in 1:(1/timestep_size)*(num_time_steps -1)){ 
constrained_pars <- c(alpha, 0, scale_parameter*(time*timestep_size)**(1/alpha), 0)

  no_hit_probability <-  no_hit_probability *(1 - stable_cdf(initial_distance + .5, constrained_pars) -   stable_cdf(initial_distance - .5, constrained_pars)*timestep_size )
            # coalescence zone has width one in our simulations

    }

constrained_probability <- 1 - no_hit_probability

#weight_list[,1] <- weight_list[,1]*constrained_probability

for (trial in 1:num_trials) {
    
   weight <- weight_list[trial, 1]
   new_filestring <- paste("trajectory_WEIGHT_alpha", alpha , "distance", initial_distance, "MU", MU, "trial", (trial -1), ".txt", sep = "") 
   sink(new_filestring)
   cat(weight)
   sink()

   #print(file[1,1])
   #print(file[2,1])
   
    
   }




#constrained_pars <- c(alpha, 0, scale_parameter*(time)**(1/alpha), 0)
extra_filestring <- paste("probability_of_constrained_trajectory_WEIGHT_alpha", alpha , "distance", initial_distance,  ".txt", sep = "") 
sink(extra_filestring)
#constrained_probability <- stable_cdf(initial_distance + .5, constrained_pars) -   stable_cdf(initial_distance - .5, constrained_pars)
cat(constrained_probability)
sink()









setwd("../..")
#system("cd ../..")