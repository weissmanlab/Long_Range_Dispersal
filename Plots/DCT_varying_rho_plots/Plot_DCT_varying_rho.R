library(libstableR)
library(argparse)
library(mdatools)
library(ggplot2)
library(gridExtra)
library(Hmisc)
library(dplyr)


for(j in 1:4){for(ZZ in 1:10){
if(j==1)
{alpha <- 1.25}
if(j==2)
{alpha <- 1.45}
if(j==3)
{alpha <- 1.65}
if(j==4)
{alpha <- 1.85}

#alpha <- 1.65
distance <- exp(ZZ)
distance_lower_bound <- .9*distance - .001
distance_upper_bound <- 1.1*distance + .001


Final_Time <-1000

Coalescence_Data_plot_ALL <-  data.frame(matrix(0, Final_Time, 3))
for( Q in 0:3)
{
rho_inverse <- 10^(-1*Q)
rho <- 1/rho_inverse



Coalescence_Data <-   read.table("ALL_runs_DCT.txt")


Coalescence_Data <- subset(Coalescence_Data, V1 == alpha)
Coalescence_Data <- subset(Coalescence_Data, V2 == 1/rho)
Coalescence_Data <- subset(Coalescence_Data, V3 > distance_lower_bound)
Coalescence_Data <- subset(Coalescence_Data, V3 < distance_upper_bound)

Coalescence_Data_Transposed <- t(Coalescence_Data)
#Coalescence_Data_Transposed <- Coalescence_Data_Transposed[-c(1:3), ]


Coalescence_Data_plot <-  data.frame(matrix(0, Final_Time, 2))
for( i in 1:Final_Time )
{Coalescence_Data_plot[i, 1] <- i -1 
 
 Coalescence_Data_plot[i, 2] <- log(Coalescence_Data_Transposed[i + 3, 1])	
 	
}

#print(Coalescence_Data_plot[, 2])
normalization_check <- sum(exp(Coalescence_Data_plot[, 2]))
#normalization_check <- sum(Coalescence_Data_plot[, 2])
#print(normalization_check)
Coalescence_Data_plot_ALL[, 1] <-Coalescence_Data_plot[, 1]
Coalescence_Data_plot_ALL[, Q+ 1] <-Coalescence_Data_plot[, 2]
}


#print(Coalescence_Data_plot)

#Coalescence_Data_plot <- subset(Coalescence_Data_plot

alpha_dummy <- 100*(alpha -1)
print(paste("Log_Plot_DCT_varying_distance_alpha_1p", alpha_dummy, "_rho_", rho, ".pdf", sep = ""))
pdf(paste("Log_Plot_DCT_varying_distance_alpha_1p", alpha_dummy, "_rho_", rho, ".pdf", sep = ""))

p <- ggplot() + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X2, color = "init dist e^01")) + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X3, color = "init dist e^02"))   + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X4, color = "init dist e^03")) + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X5, color = "init dist e^04"))   +  geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X6, color = "init dist e^05")) + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X7, color = "init dist e^06")) + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X8, color = "init dist e^07")) + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X9, color = "init dist e^08")) +  geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X10, color = "init dist e^09")) + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X11, color = "init dist e^10")) + geom_point() + geom_point() + labs( x = "Time", y ="Log Dist of Coalescence Times") + ggtitle(paste("Alpha", alpha, "Rho", 1/rho_inverse))#+ labs( x = "Time", y ="Log Dist of Coalescence Times") + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X12, color = "init dist 11")) 


print(p)

dev.off()

}}