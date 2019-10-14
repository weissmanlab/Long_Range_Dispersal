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
#CDF_table <-data.frame(read.table("10_6_runs_alpha_p75_CDF_of_coalescence_timesalpha_value_0.75distance_value_0000000rho_inverse_1e-05.txt"))
CDF_table <-data.frame(read.table("million_step_high_alpha_c_5_CDF.txt"))
scale_parameter <- 5
alpha <- 2.05
rho <- 1
#rho <- 100000
rho_inverse <- 1/rho
CDF_table <- subset(CDF_table, V1 == alpha)
if(rho ==1){CDF_table <- subset(CDF_table, V2 == rho_inverse )}
if(rho_inverse ==.00001){CDF_table <- subset(CDF_table, V2 == "1e-05" )}
#alpha <- .50
CDF_table <- CDF_table[-14, ]
#CDF_table <- CDF_table[-1, ]

#CDF_table <- CDF_table[-1, ]

#print(Log_Comp_CDF_table[,1])

D <- ((scale_parameter)^(alpha))/2
improper_norm <- 1/(2*pi*rho*(1- alpha)*D + 1)
improper_norm <- max(CDF_table[,5])
#print(improper_norm)
#print(max(CDF_table[,5]))
#print(log(improper_norm))
if(alpha < 1)
{
CDF_table[,3] <- log10(CDF_table[,3] )

#shifted CDF
#CDF_table[,5] <- CDF_table[,5] - improper_norm 
#CDF_table[,6] <-  CDF_table[,6] - improper_norm 
#CDF_table[,7] <-  CDF_table[,7] - improper_norm 
#CDF_table <- CDF_table[-1, ]

CDF_table[,5] <- log10( CDF_table[,5] )
CDF_table[,6] <- log10( CDF_table[,6] )
CDF_table[,7] <- log10(CDF_table[,7] )

fun.1 <- function(x) log10(improper_norm - 10^((1- 1/alpha)*x + log10(gamma(1+1/alpha)/(((2*D)^(1/alpha))*pi*rho))  + log10(1/(1/alpha - 1))) ) 
p <- ggplot(CDF_table, aes(V3, V5)) +
  geom_point() +  geom_pointrange(data=CDF_table, aes(x = V3, y =V5, ymin =V7, ymax =V6), pch  = 0) +
stat_function(fun=fun.1) + labs( x = "Log Time", y ="Log CDF") + ggtitle(paste("Alpha", alpha, "Distance 0")) + xlim(0, 6.2)

  
  print(p)
}




if(alpha >1)
{
CDF_table[,3] <- log10(CDF_table[,3] )
CDF_table[,5] <- log10(CDF_table[,5] )
CDF_table[,6] <- log10(CDF_table[,6] )
CDF_table[,7] <- log10(CDF_table[,7] )

fun.1 <- function(x)   log10(1 - 10^((1/alpha -1)*x + log10((alpha)*rho*((2*D)^(1/alpha))*sin(pi/alpha)))) #log((alpha -1)*rho*(2*D)^(1/alpha)*sin(pi/alpha)) + log(1/(1-1/alpha))  
#if( alpha >= 2){fun.1 <- function(x)   (.5)*x }#+ log10((alpha)*rho*((2*D)^(1/alpha))*sin(pi/alpha))}
p <- ggplot(CDF_table, aes(V3, V5)) +
  geom_point() +  geom_pointrange(data=CDF_table, aes(x = V3, y =V5, ymin =V7, ymax =V6), pch  = 0) +
stat_function(fun=fun.1) + labs( x = "Log Time", y ="Log CDF") + ggtitle(paste("Alpha", alpha, "Distance 0"))

 
  
  print(p)
}




if(alpha == 1)
{
CDF_table[,3] <- log(CDF_table[,3] )
CDF_table[,5] <- log(CDF_table[,5] )
CDF_table[,6] <- log(CDF_table[,6] )
CDF_table[,7] <- log(CDF_table[,7] )


fun.1 <- function(x)  log(2*pi*rho*D) -log(log(10)*x + log(2*D))    #log(2*pi*rho*D/( 2*pi*rho*D+ log(1 +D*2*exp(x))))
p <- ggplot(CDF_table, aes(V3, V5)) +
  geom_point() +  geom_pointrange(data=CDF_table, aes(x = V3, y =V5, ymin =V7, ymax =V6), pch  = 0) +
stat_function(fun=fun.1) + labs( x = "Log Time", y ="Log CDF") + ggtitle(paste("Alpha", alpha, "Distance 0"))
 
  
  print(p)
}