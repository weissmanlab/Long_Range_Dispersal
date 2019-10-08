library(libstableR)
library(argparse)
library(mdatools)
library(ggplot2)
library(gridExtra)
library(Hmisc)
library(dplyr)
library(ggpmisc)
#Final_Time <- 10000
CDF_table <-data.frame(read.table("CDF_of_coalescence_timesalpha_value_0.23distance_value_0000000rho_inverse_1e-05.txt"))
CDF_table <- CDF_table[-1, ]
alpha <- .23
#print(Log_Comp_CDF_table[,1])
scale_parameter <- 1
D <- ((scale_parameter)^(1/alpha))/2
rho <- 100000
improper_norm <- 1/(2*pi*rho*(1- alpha)*D + 1)

CDF_table[,3] <- log(CDF_table[,3] )
CDF_table[,5] <- log(improper_norm - CDF_table[,5] )
CDF_table[,6] <- log(improper_norm - CDF_table[,6] )
CDF_table[,7] <- log(improper_norm - CDF_table[,7] )

p <- ggplot(CDF_table, aes(V3, V5)) +
  geom_point() + geom_abline(intercept = -11.3, slope = (1-1/alpha), color="red", 
                 linetype="dashed", size=1.5) + geom_pointrange(data=CDF_table, aes(x = V3, y =V5, ymin =V7, ymax =V6), pch  = 0)
 
  
  print(p)
