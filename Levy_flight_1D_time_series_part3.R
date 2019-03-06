library(stabledist)
ALPHA = 1.25
rho_inverse = 1
distance = 0
scale_parameter = 250




#string = paste("./rho_inverse_", rho_inverse, "alpha_", alpha, sep = "")
string = paste("./alpha_value_", ALPHA, "/distance_value_", distance, sep = "")
#print(string)
#setwd(string)


#print(Sys.glob(paste("histogram_of_boot*", rho_inverse, "*", sep = "")))

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
Bootstrapped_Mean_Homozygosities_Histogram <-  read.table(Sys.glob(paste("histogram_of_boot*", rho_inverse, "*", sep = "")))
Single_Trial_Homozygosities_Histogram <-   read.table(Sys.glob(paste("histogram_of_single*", rho_inverse, "*", sep = "")))
}
if(rho_inverse == 1)
{
Bootstrapped_Mean_Homozygosities_Histogram <-  read.table(Sys.glob(paste("histogram_of_boot*", "*rho_inverse_1.txt", sep = "")))
Single_Trial_Homozygosities_Histogram <-   read.table(Sys.glob(paste("histogram_of_single*",  "*rho_inverse_1.txt", sep = "")))
}


	Bootstrapped_Mean_Homozygosities_Histogram <- subset(Bootstrapped_Mean_Homozygosities_Histogram, V2 == mu)
	
	Single_Trial_Homozygosities_Histogram  <- subset(Single_Trial_Homozygosities_Histogram , V2 == mu)
 #print(sum(Single_Trial_Homozygosities_Histogram[,4]))
 #print(sum(Bootstrapped_Mean_Homozygosities_Histogram[,4]))


dummy_sum <- 0
for(N in 1:100)
{

dummy_prob <- exp(-mu*N)*dstable(x = distance, alpha = ALPHA , beta = 0, gamma = scale_parameter*(N**(1/ALPHA)), delta = 0, pm = 0,log = FALSE,tol = 64*.Machine$double.eps, zeta.tol = NULL,subdivisions = 1000)
dummy_value <- Numerical_MH_distance_zero[, 5]*exp(-mu*N)
#resulting_homozygosity <- 
#make sure this scale parameter is right - no extra factor of 2

#dummy2 <- dstable(x, alpha, beta, gamma)
dummy_sum <- dummy_sum + dummy_prob

}


Original_Histogram <- subset(Bootstrapped_Mean_Homozygosities_Histogram, select = c(V3, V4) )

#Original_Histogram_Mean <- 0


print(Smoothed_Histogram)




for(N in 1:100)
{

dummy_prob <- exp(-mu*N)*dstable(x = distance, alpha = ALPHA , beta = 0, gamma = scale_parameter*(N**(1/ALPHA)), delta = 0, pm = 0,log = FALSE,tol = 64*.Machine$double.eps, zeta.tol = NULL,subdivisions = 1000)
dummy_value <- Numerical_MH_distance_zero[, 5]*exp(-mu*N)
#resulting_homozygosity <- 
#make sure this scale parameter is right - no extra factor of 2


}




#print(dummy_sum)
setwd("./../../")

}






