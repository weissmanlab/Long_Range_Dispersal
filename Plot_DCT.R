library(libstableR)
library(argparse)
library(mdatools)
library(ggplot2)
library(gridExtra)
library(Hmisc)
library(dplyr)

alpha <- 1.65
rho_inverse <- .1
rho <- 1/rho_inverse

Final_Time <-1000

Coalescence_Data_plot_ALL <-  data.frame(matrix(0, Final_Time, 11))
for( Q in 1:10)
{
distance <- exp(Q)
distance_lower_bound <- .9*distance - .001
distance_upper_bound <- 1.1*distance + .001

Coalescence_Data <-   read.table("ALL_runs_DCT.txt")


Coalescence_Data <- subset(Coalescence_Data, V1 == alpha)
Coalescence_Data <- subset(Coalescence_Data, V2 == rho_inverse)
Coalescence_Data <- subset(Coalescence_Data, V3 > distance_lower_bound)
Coalescence_Data <- subset(Coalescence_Data, V3 < distance_upper_bound)

Coalescence_Data_Transposed <- t(Coalescence_Data)
#Coalescence_Data_Transposed <- Coalescence_Data_Transposed[-c(1:3), ]


Coalescence_Data_plot <-  data.frame(matrix(0, Final_Time, 2))
for( i in 1:Final_Time )
{Coalescence_Data_plot[i, 1] <- i -1 
 
 Coalescence_Data_plot[i, 2] <- log(Coalescence_Data_Transposed[i + 3, 1])	
 	
}

print(Coalescence_Data_plot[, 2])
normalization_check <- sum(exp(Coalescence_Data_plot[, 2]))
#normalization_check <- sum(Coalescence_Data_plot[, 2])
#print(normalization_check)
Coalescence_Data_plot_ALL[, 1] <-Coalescence_Data_plot[, 1]
Coalescence_Data_plot_ALL[, Q+ 1] <-Coalescence_Data_plot[, 2]
}


#print(Coalescence_Data_plot)

#Coalescence_Data_plot <- subset(Coalescence_Data_plot

p <- ggplot() + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X2)) + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X6, color = "init dist 1")) + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X7, color = "init dist 2")) + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X8, color = "init dist 3")) + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X9, color = "init dist 4"))   + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X10, color = "init dist 5")) + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X11, color = "init dist 6")) + geom_point()

print(p)