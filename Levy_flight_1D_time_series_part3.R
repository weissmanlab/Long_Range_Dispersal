library(stabledist)
library(data.table)




ALPHA = 1.25
rho_inverse = 1
distance = 0
scale_parameter = 250
string = paste("./alpha_value_", ALPHA, "/distance_value_", distance, sep = "")
string2 = paste("smoothed_mean_homozygosity_alpha_value_", ALPHA, "distance_value_", distance, "rho_inverse_",rho_inverse, ".txt", sep = "")
sink(string2)




for(k in 1:3)
{

   if(k == 1) {mu = 0.001}
   if(k == 2) {mu = 0.01}
   if(k == 3) {mu = 0.1}
   
  Numerical_MH_distance_zero <- read.table("DF_ALL_zeros_numerics.txt")

Numerical_MH_distance_zero <- subset(Numerical_MH_distance_zero, V1 == ALPHA)
Numerical_MH_distance_zero <- subset(Numerical_MH_distance_zero, V2 == rho_inverse)
Numerical_MH_distance_zero <- subset(Numerical_MH_distance_zero, V3 == mu)
#print(Numerical_MH_distance_zero_sub3[,5])
setwd(string)

if(rho_inverse != 1)
{
#Bootstrapped_Mean_Homozygosities_Histogram <-  read.table(Sys.glob(paste("histogram_of_boot*", rho_inverse, "*", sep = "")))
#Single_Trial_Homozygosities_Histogram <-   read.table(Sys.glob(paste("histogram_of_single*", rho_inverse, "*", sep = "")))
Sorted_List_of_Bootstrapped_Mean_Homozygosities <- read.table(Sys.glob(paste("sorted_list_of_boot*", rho_inverse, "*", sep = "")))
## sorted list has weights in in its last columns, so it can serve as a continuous histogram without rounding/binning error.
}
if(rho_inverse == 1)
{
#Bootstrapped_Mean_Homozygosities_Histogram <-  read.table(Sys.glob(paste("histogram_of_boot*", "*rho_inverse_1.txt", sep = "")))
#Single_Trial_Homozygosities_Histogram <-   read.table(Sys.glob(paste("histogram_of_single*",  "*rho_inverse_1.txt", sep = "")))
Sorted_List_of_Bootstrapped_Mean_Homozygosities <- read.table(Sys.glob(paste("sorted_list_of_boot*", "*rho_inverse_1.txt", sep = "")))
## sorted list has weights in in its last columns, so it can serve as a continuous histogram without rounding/binning error.
}


	Sorted_List_of_Bootstrapped_Mean_Homozygosities <- subset(Sorted_List_of_Bootstrapped_Mean_Homozygosities, V2 == mu)
	
	
	Original_Mean <- sum(Sorted_List_of_Bootstrapped_Mean_Homozygosities[,3]*Sorted_List_of_Bootstrapped_Mean_Homozygosities[,4]) 
	
	
	#Sorted_List_of_Bootstrapped_Mean_Homozygosities  <- subset(Sorted_List_of_Bootstrapped_Mean_Homozygosities , V2 == mu)
 #print(sum(Single_Trial_Homozygosities_Histogram[,4]))
 #print(sum(Bootstrapped_Mean_Homozygosities_Histogram[,4]))


dummy_sum <- 0
for(N in 1:100)
{

dummy_prob <- exp(-mu*N)*dstable(x = distance, alpha = ALPHA , beta = 0, gamma = scale_parameter*(N**(1/ALPHA)), delta = 0, pm = 0,log = FALSE,tol = 64*.Machine$double.eps, zeta.tol = NULL,subdivisions = 1000)
dummy_value <- Numerical_MH_distance_zero[, 5]*exp(-mu*N)

#make sure this scale parameter is right - no extra factor of 2

dummy_sum <- dummy_sum + dummy_prob

}


Original_Sorted_List <- subset(Sorted_List_of_Bootstrapped_Mean_Homozygosities, select = c(V3, V4) )
Smoothed_Sorted_List <- (1-dummy_sum)*Original_Sorted_List
setnames(Smoothed_Sorted_List, old = c('V3', 'V4'),  new = c('sampled_MH_value', 'probability'))

#We rescale the entries of the original list so that the new sum will be normalized.  The rescaling factor (1-dummy_sum) is very close to 1, around .99


#Original_Histogram <- subset(Sorted_List_of_Bootstrapped_Mean_Homozygosities, select = c(V3, V4) )


#Original_Histogram_Mean <- 0







for(N in 1:100)
{

probability <- exp(-mu*N)*dstable(x = distance, alpha = ALPHA , beta = 0, gamma = scale_parameter*(N**(1/ALPHA)), delta = 0, pm = 0,log = FALSE,tol = 64*.Machine$double.eps, zeta.tol = NULL,subdivisions = 1000)
sampled_MH_value <- Numerical_MH_distance_zero[, 5]*exp(-mu*N)

#make sure this scale parameter is right - no extra factor of 2
df <- data.frame( sampled_MH_value, probability)
if(N ==1)
{Extra_Sorted_List <- df} 
if(N !=1)
{Extra_Sorted_List <- rbind(Extra_Sorted_List, df)} 

}

#print(Extra_Sorted_List)

# Now just combine and bubble sort the two lists - No need for merge sort since they're relatively short
Total_List = rbind(Smoothed_Sorted_List, Extra_Sorted_List)


Smoothed_Sorted_List <- Total_List[with(Total_List, order(sampled_MH_value)),]


CI_lower_sum <-0
CI_upper_sum <-0

CI_lower_index <-0
CI_upper_index <-0

CI_lower_BOOL <-0
CI_upper_BOOL <-0



for (N in 1:nrow(Smoothed_Sorted_List))
{  if( CI_lower_sum < .2)
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
#print(Smoothed_Sorted_List[CI_lower_index,])
#print(Smoothed_Sorted_List[CI_upper_index,])


print(paste(distance, mu, Original_Mean, Smoothed_Sorted_List[CI_lower_index,1], Smoothed_Sorted_List[CI_upper_index,1]))

#convert pdf to cdf


#print(Smoothed_Sorted_List[,2])


#print(sum(Smoothed_Sorted_List[,1]*Smoothed_Sorted_List[,2]))


#print(dummy_sum)
setwd("./../../")

}

sink()




