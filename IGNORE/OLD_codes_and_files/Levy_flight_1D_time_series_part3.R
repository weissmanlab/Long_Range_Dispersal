library(stabledist)
library(data.table)

ALPHA = 1.25
rho_inverse = 1
distance = 5.0
#negative_distance = -1.0*distance
scale_parameter = 250
string = paste("./alpha_value_", ALPHA, "/distance_value_", distance, sep = "")
string2 = paste("smoothed_mean_homozygosity_alpha_value_", ALPHA, "distance_value_", distance, "rho_inverse_",rho_inverse, ".txt", sep = "")
sink(string2)


for(k in 1:3)
{  if(k == 1) {mu = 0.001}
   if(k == 2) {mu = 0.01}
   if(k == 3) {mu = 0.1}
   
  Numerical_MH_distance_zero <- read.table("DF_ALL_zeros_numerics.txt")

Numerical_MH_distance_zero <- subset(Numerical_MH_distance_zero, V1 == ALPHA)
Numerical_MH_distance_zero <- subset(Numerical_MH_distance_zero, V2 == rho_inverse)
Numerical_MH_distance_zero <- subset(Numerical_MH_distance_zero, V3 == mu)
#print(Numerical_MH_distance_zero_sub3[,5])
setwd(string)

# set working directory

if(rho_inverse != 1)
{Sorted_List_of_Bootstrapped_Mean_Homozygosities <- read.table(Sys.glob(paste("sorted_list_of_boot*", rho_inverse, "*", sep = "")))
}
if(rho_inverse == 1)
{Sorted_List_of_Bootstrapped_Mean_Homozygosities <- read.table(Sys.glob(paste("sorted_list_of_boot*", "*rho_inverse_1.txt", sep = "")))
}
Sorted_List_of_Bootstrapped_Mean_Homozygosities <- subset(Sorted_List_of_Bootstrapped_Mean_Homozygosities, V2 == mu)
Original_Mean <- sum(Sorted_List_of_Bootstrapped_Mean_Homozygosities[,3]*Sorted_List_of_Bootstrapped_Mean_Homozygosities[,4]) 
	
	## We read in the sorted list of bootstrapped mean homozygosities and compute the overall mean

############################### Now we smooth the data


q_small <- -1.0*(distance/10.0):1.0*(distance/10.0) # This is the length scale very close to the origin that trajectories rarely hit and is thus undersampled
q_large <- -1.0*distance:1.0*distance # This is the length scale close enough to the origin to serve as a cutoff.  Paths that have crossed this bound are unlikely to hit origin on relevant timescales
prob_sum <- 0


for(N in 1:as.integer(10.0/mu))
{dummy_prob <- pstable(q_small, alpha = ALPHA , beta = 0, gamma = scale_parameter, delta = 0, pm = 0,log = FALSE,tol = 64*.Machine$double.eps,subdivisions = 1000)*pstable(q_large, alpha = ALPHA , beta = 0, gamma = scale_parameter, delta = 0, pm = 0,log = FALSE,tol = 64*.Machine$double.eps, subdivisions = 1000)**(N-1)
# we apporximate the probability of hitting (near) the origin at a given timestep as the quantile from -qsmall to qsmall times the probability that the path hasn't strayed extremely far from the origin at a previous timestep, which we quantify as the quantile from -qlarge to qlarge to the N-1 power.
dummy_value <- Numerical_MH_distance_zero[, 5]*exp(-2*mu*N)
prob_sum <- prob_sum + dummy_prob
}
## In the above loop we determine the probability of any rare event that was likely undersampled - this is quantified by prob_sum
Original_Sorted_List <- subset(Sorted_List_of_Bootstrapped_Mean_Homozygosities, select = c(V3, V4) )
Smoothed_Sorted_List <- Original_Sorted_List
Smoothed_Sorted_List[, 1] <-(1 - prob_sum)*Smoothed_Sorted_List[, 1] 
#setnames(Smoothed_Sorted_List, old = c('V3', 'V4'),  new = c('sampled_MH_value', 'probability'))

#We rescale the entries of the original list so that the new sum will be normalized.  

# Now we can add the sampled events to the rare events that were likely undersampled in an effort to smooth every entry in the sorted list.  

for(N in 1:as.integer(10.0/mu))
{probability <-  pstable(q_small, alpha = ALPHA , beta = 0, gamma = scale_parameter, delta = 0, pm = 0,log = FALSE,tol = 64*.Machine$double.eps,subdivisions = 1000)*pstable(q_large, alpha = ALPHA , beta = 0, gamma = scale_parameter, delta = 0, pm = 0,log = FALSE,tol = 64*.Machine$double.eps,subdivisions = 1000)**(N-1)
# we apporximate the probability of hitting (near) the origin at a given timestep as the quantile from -qsmall to qsmall times the probability that the path hasn't strayed extremely far from the origin at a previous timestep, which we quantify as the quantile from -qlarge to qlarge to the N-1 power.
sampled_MH_value <- Numerical_MH_distance_zero[, 5]*exp(-2*mu*N)
Smoothed_Sorted_List[, 1] <- Smoothed_Sorted_List[, 1] + probability*sampled_MH_value
}

## The above loop smoothes all the entries in the list of bootstrapped mean homozygosities by introducing the effects of rare but important paths that likely were not sampled.
#We smooth every bootstrapped mean homozygosity to account for rare and undersampled trajectories.

############################### Above we smooth the data


####################################### Below we compute the confidence interva;
CI_lower_sum <-0
CI_upper_sum <-0

CI_lower_index <-0
CI_upper_index <-0

CI_lower_BOOL <-0.  ####### Confidence intervals for new smoothed data
CI_upper_BOOL <-0


CI_lower_sum_OLD <-0
CI_upper_sum_OLD <-0

CI_lower_index_OLD <-0
CI_upper_index_OLD <-0.   ######### Confidence intervals for old, unsmoothed data

CI_lower_BOOL_OLD <-0
CI_upper_BOOL_OLD <-0

######### Compute old and new confidence intervals
for (N in 1:nrow(Smoothed_Sorted_List))
{  
	if( CI_lower_sum < .2)
	{CI_lower_sum <- CI_lower_sum + Smoothed_Sorted_List[N, 2]}
   if(CI_lower_sum >= .2 & CI_lower_BOOL ==0)
   {CI_lower_index <- N
   	CI_lower_BOOL <- 1  }
	
	
	if( CI_upper_sum < .8)
	{CI_upper_sum <- CI_upper_sum + Smoothed_Sorted_List[N, 2]}
	if(CI_upper_sum >= .8 & CI_upper_BOOL ==0)
   {CI_upper_index <- N
   	CI_upper_BOOL <- 1  }

}



for (N in 1:nrow(Sorted_List_of_Bootstrapped_Mean_Homozygosities))
{  if( CI_lower_sum_OLD < .2)
	{CI_lower_sum_OLD <- CI_lower_sum_OLD + Original_Sorted_List[N, 2]}
   if(CI_lower_sum_OLD >= .2 & CI_lower_BOOL_OLD ==0)
   {CI_lower_index_OLD <- N
   	CI_lower_BOOL_OLD <- 1  }
	
	
	if( CI_upper_sum_OLD < .8)
	{CI_upper_sum_OLD <- CI_upper_sum_OLD + Original_Sorted_List[N, 2]}
	if(CI_upper_sum_OLD >= .8 & CI_upper_BOOL_OLD ==0)
   {CI_upper_index_OLD <- N
   	CI_upper_BOOL_OLD <- 1  }

}

##### True confidence intervals will be the lower of the smoothed and unsmoothed lower bounds and the higher of the smoothed and unsmoothed upper bounds

True_CI_lower <- 0

if(CI_lower_index_OLD <=  CI_lower_index)
{True_CI_lower <- Original_Sorted_List[CI_lower_index_OLD, 1]  }

if(CI_lower_index_OLD >  CI_lower_index)
{True_CI_lower <- Smoothed_Sorted_List[CI_lower_index, 1]  }


True_CI_upper <- 0

if(CI_upper_index_OLD >=  CI_upper_index)
{True_CI_upper <- Original_Sorted_List[CI_upper_index_OLD, 1]  }

if(CI_upper_index_OLD <  CI_upper_index)
{True_CI_upper <- Smoothed_Sorted_List[CI_upper_index, 1]  }


cat(noquote(paste(distance, mu, Original_Mean, True_CI_lower, True_CI_upper)))
cat("\n")

#print(dummy_sum)
setwd("./../../")

}

sink()




