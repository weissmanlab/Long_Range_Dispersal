#########################################################
#install.packages("argparse")
library(libstableR)
library(argparse)
library(mdatools)
args <- commandArgs(trailingOnly = TRUE)
if (length(args)!=5) {

  print(length(args))
  stop("Five input arguments are required - alpha, initial distance, number of trials for each time, max number of timesteps, scale parameter", call.=FALSE)
  
}


alpha <- as.double(args[1])
initial_distance <- as.double(args[2])
num_trials <- as.integer(args[3]) # number of trials for each time
max_num_time_steps <- as.integer(args[4]) # number of times considered
scale_parameter <- as.double(args[5])
pars <- c(alpha, 0, scale_parameter, 0)
# store name of directory and change into said directory
dir <- paste("./alpha_value_", alpha, "/distance_value_", initial_distance, sep = "")
#print(dir)
setwd(dir)
#system(cd_to_dir)

for (time in 1:(max_num_time_steps -1)) {for (trial in 1:num_trials) {
    
   filestring <- paste("free_and_contrained_CHOSEN_TIME_jump_sizesalpha", alpha , "distance", initial_distance, "num_time_steps", time , "trial", (trial -1), ".txt", sep = "") 
   #print(filestring)
   # -1 adjusts for fact that c++ loops from part 1 start at zero and r loops start at 1
   file <- read.table(filestring)
   free_step <- file[1,1]
   constrained_step <- file[2,1]
   
   weight <- stable_pdf(constrained_step, pars)/stable_pdf(free_step , pars)
   
   new_filestring <- paste("trajectory_WEIGHT_alpha", alpha , "distance", initial_distance, "num_time_steps", time , "trial", (trial -1), ".txt", sep = "") 
   sink(new_filestring)
   print(weight)
   sink()

   #print(file[1,1])
   #print(file[2,1])
   
    
   }
constrained_pars <- c(alpha, 0, scale_parameter*(time)**(1/alpha), 0)
extra_filestring <- paste("probability_of_constrained_trajectory_WEIGHT_alpha", alpha , "distance", initial_distance, "num_time_steps", time , ".txt", sep = "") 
sink(extra_filestring)
constrained_probability <- stable_cdf(initial_distance + .5, constrained_pars) -   stable_cdf(initial_distance - .5, constrained_pars)
cat(constrained_probability)
sink()
}





setwd("../..")
#system("cd ../..")