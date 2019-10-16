library(libstableR)
library(argparse)
library(mdatools)
library(ggplot2)
library(gridExtra)
library(Hmisc)
library(dplyr)
library(ggpmisc)
#Final_Time <- 10000
#CDF_table <-data.frame(read.table("CDF_dist0_c1.txt"))
#CDF_table <-data.frame(read.table("million_step_high_alpha_c_10_CDF.txt"))
#CDF_table <-data.frame(read.table("million_step_high_alpha_c_5_CDF.txt"))
#CDF_table <-data.frame(read.table("million_step_alpha_1_c_1_CDF.txt"))
#CDF_table <-data.frame(read.table("10_6_runs_alpha_p75_CDF_of_coalescence_timesalpha_value_0.75distance_value_0000000rho_inverse_1e-05.txt"))
#CDF_table <-data.frame(read.table("CDF_of_coalescence_timesalpha_value_0.5distance_value_0000000rho_inverse_1e-05.txt"))
#CDF_table <-data.frame(read.table("CDF_of_coalescence_timesalpha_value_0.75distance_value_0000000rho_inverse_1e-05.txt"))
#CDF_table <-data.frame(read.table("CDF_of_coalescence_timesalpha_value_1distance_value_0000000rho_inverse_1000.txt"))
#CDF_table <-data.frame(read.table("CDF_runs_alpha_1_dist0_c_250.txt"))
#CDF_table <-data.frame(read.table("CDF_runs_high_alpha_c_5.txt"))
CDF_table <-data.frame(read.table("CDF_runs_alpha_1_c_1.txt"))
#CDF_table <-data.frame(read.table("CDF_runs_low_alpha_c1.txt"))
scale_parameter <- 1
alpha <- 1
#rho <- 100000
rho <- 10
#rho <- 1
rho_inverse <- 1/rho
CDF_table <- subset(CDF_table, V1 == alpha)
CDF_table <- subset(CDF_table, V2 == rho_inverse )
#if(rho_inverse ==.00001){CDF_table <- subset(CDF_table, V2 == "1e-05" )}
#if(rho_inverse ==.00001){CDF_table[,2] <- .00001}
#alpha <- 1.0
#CDF_table <- CDF_table[-14, ]
CDF_table <- CDF_table[-1, ]

#CDF_table <- CDF_table[-1, ]

#print(Log_Comp_CDF_table[,1])
#alpha <- 2.0

D <- ((scale_parameter)^(alpha))/2
#improper_norm <- 1/(2*pi*rho*(1- alpha)*D + 1)
#improper_norm <- improper_norm/.985
#print(improper_norm)
improper_norm <- max(CDF_table[,5]) #/.999
print(improper_norm)
#print(max(CDF_table[,5])/.99)

if(alpha < 1)
{
CDF_table[,3] <- CDF_table[,3]*(2*((1/CDF_table[,2])^alpha)*D)^(1/(1 - alpha))

CDF_table[,6] <-  max(CDF_table[,5]) - CDF_table[,6] 
CDF_table[,7] <-  max(CDF_table[,5]) - CDF_table[,7] 
CDF_table[,5] <-  max(CDF_table[,5])- CDF_table[,5] 



CDF_table[,5] <- CDF_table[,5]*pi*(1/alpha -1)/gamma(1/alpha +1)
CDF_table[,6] <- CDF_table[,6]*pi*(1/alpha -1)/gamma(1/alpha +1)
CDF_table[,7] <- CDF_table[,7]*pi*(1/alpha -1)/gamma(1/alpha +1)



CDF_table[,3] <- log10(CDF_table[,3] )
CDF_table[,5] <- log10(CDF_table[,5] )
CDF_table[,6] <- log10(CDF_table[,6] )
CDF_table[,7] <- log10(CDF_table[,7] )
fun.1 <- function(x) (1- 1/alpha)*x  
p <- ggplot(CDF_table, aes(V3, V5)) +
  geom_point() +  geom_pointrange(data=CDF_table, aes(x = V3, y =V5, ymin =V7, ymax =V6), pch  = 0) +
stat_function(fun=fun.1) + labs( x = "Log Time", y ="Log Comp CDF") + ggtitle(paste("Alpha", alpha, "Distance 0")) 
  
  print(p)
}




if(alpha >1)
{

CDF_table[,3] <- CDF_table[,3]*(2*((1/CDF_table[,2])^alpha)*D)^(1/(1 - alpha))


CDF_table[,5] <- 1 - CDF_table[,5] 
CDF_table[,6] <- 1 - CDF_table[,6] 
CDF_table[,7] <- 1 - CDF_table[,7] 



CDF_table[,5] <- CDF_table[,5]/(alpha*sin(pi/alpha))
CDF_table[,6] <- CDF_table[,6]/(alpha*sin(pi/alpha))
CDF_table[,7] <- CDF_table[,7]/(alpha*sin(pi/alpha))



CDF_table[,3] <- log10(CDF_table[,3] )
CDF_table[,5] <- log10(CDF_table[,5] )
CDF_table[,6] <- log10(CDF_table[,6] )
CDF_table[,7] <- log10(CDF_table[,7] )

fun.1 <- function(x)   ( 1/alpha - 1)*x 

p <- ggplot(CDF_table, aes(V3, V5)) +
  geom_point() +  geom_pointrange(data=CDF_table, aes(x = V3, y =V5, ymin =V7, ymax =V6), pch  = 0) +
stat_function(fun=fun.1) + labs( x = "Log Scaled Time", y ="Log Scaled Comp CDF") + ggtitle(paste("Alpha", alpha, "Distance 0"))

 
  
  print(p)
}




if(alpha == 1)
{

CDF_table[,3] <- CDF_table[,3]*(2*D)*exp(2*pi*rho*D)

CDF_table[,5] <- 1 - CDF_table[,5] 
CDF_table[,6] <- 1 - CDF_table[,6] 
CDF_table[,7] <- 1 - CDF_table[,7] 



CDF_table[,5] <- CDF_table[,5]/(2*pi*rho*D)
CDF_table[,6] <- CDF_table[,6]/(2*pi*rho*D)
CDF_table[,7] <- CDF_table[,7]/(2*pi*rho*D)



CDF_table[,3] <- log10(CDF_table[,3] )
CDF_table[,5] <- log10(CDF_table[,5] )
CDF_table[,6] <- log10(CDF_table[,6] )
CDF_table[,7] <- log10(CDF_table[,7] )



#fun.1 <- function(x)      log10(2*pi*rho*D) -log10(log(10)*x + log(2*D)  + 2*pi*rho*D ) #log10(2*pi*rho*D) -log10(log(10)*x + log(10)*log(2*D) + 2*pi*rho*D )    
fun.1 <- function(x)      - log10(log(10)*x) 
#log(2*pi*rho*D/( 2*pi*rho*D+ log(1 +D*2*exp(x))))
p <- ggplot(CDF_table, aes(V3, V5)) +
  geom_point() +  geom_pointrange(data=CDF_table, aes(x = V3, y =V5, ymin =V7, ymax =V6), pch  = 0) +
stat_function(fun=fun.1) + labs( x = "Log Time", y ="Log Comp CDF") + ggtitle(paste("Alpha", alpha, "Distance 0"))
 
  
  print(p)
}




