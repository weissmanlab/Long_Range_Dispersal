library(libstableR)
library(argparse)
library(mdatools)
library(ggplot2)
library(gridExtra)
library(Hmisc)
library(dplyr)
library(ggpmisc)
#Final_Time <- 10000
Log_Comp_CDF_table <-data.frame(read.table("Log_comp_CDF_of_coalescence_timesalpha_value_0.83distance_value_0000000rho_inverse_1.txt"))
Log_Comp_CDF_table <- Log_Comp_CDF_table[-1, ]
alpha <- .83
print(Log_Comp_CDF_table[,1])


p <- ggplot(Log_Comp_CDF_table, aes(V3, V5)) +
  geom_point() + geom_abline(intercept = .00, slope = (1-1/alpha), color="red", 
                 linetype="dashed", size=1.5) + geom_pointrange(data=Log_Comp_CDF_table, aes(x = V3, y =V5, ymin =V6, ymax =V7), pch  = 0)
 
  
  print(p)
