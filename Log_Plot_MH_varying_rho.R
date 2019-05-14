library(ggplot2)
library(gridExtra)
library(Hmisc)
library(dplyr)
alpha <- 1.65
rho_inverse <- .1
mu <- .1
#par(mfrow=c(3,3))
x <- list()
y <- list()
p <- list()
k <- 1

for(j in 1:4){#for(i in 0:2){
if(j==1)
{alpha <- 1.25}
if(j==2)
{alpha <- 1.45}
if(j==3)
{alpha <- 1.65}
if(j==4)
{alpha <- 1.85}
Simulation_Data <-   read.table("ALL_runs_MH.txt")
Semianalytic_Data <-read.table("numeric_MH_LOG_SCALE_ALL_but_zeros.txt")

Simulation_Data <- Simulation_Data[!(Simulation_Data$V3 == 1 ),]
#convert simulation data to log scale
Simulation_Data[,4] <-log(Simulation_Data[,4])
Simulation_Data[,5] <-log(Simulation_Data[,5])
Simulation_Data[,6] <-log(Simulation_Data[,6])
Simulation_Data[,7] <-log(Simulation_Data[,7])

Simulation_Data <- subset(Simulation_Data, V3 < 1)
Semianalytic_Data  <- subset(Semianalytic_Data, V3 <1)
#Simulation_Data <- Simulation_Data[!(Simulation_Data$V7 - Simulation_Data$V6  > 5),]
#Semianalytic_Data <- Semianalytic_Data[!(Semianalytic_Data$V7 - Semianalytic_Data$V6  > 5),]

#Simulation_Data <- Simulation_Data[!(Simulation_Data$V5   < -30),]

#Simulation_Data  <- subset(Simulation_Data, V1 < 1.8)





Simulation_Data  <- subset(Simulation_Data, V1 == alpha)
#Simulation_Data  <-subset(Simulation_Data, V2 ==rho_inverse)
Simulation_Data <- subset(Simulation_Data, V3 ==mu)
Simulation_Data <- subset(Simulation_Data, V4 < 10.8)



#Semianalytic_Data  <- subset(Semianalytic_Data, V1 < 1.8)

Semianalytic_Data  <- subset(Semianalytic_Data, V1 == alpha)
#Semianalytic_Data  <-subset(Semianalytic_Data, V2 ==rho_inverse)
Semianalytic_Data  <- subset(Semianalytic_Data, V3 ==mu)
Semianalytic_Data <- subset(Semianalytic_Data, V4 < 10.8)
if(alpha == 1.65){Simulation_Data <- subset(Simulation_Data, V4 < 9.8)}
if(alpha == 1.65){Semianalytic_Data <- subset(Semianalytic_Data, V4 < 9.8)}

if(alpha == 1.85){Simulation_Data <- subset(Simulation_Data, V4 < 8.8)}
if(alpha == 1.85){Semianalytic_Data <- subset(Semianalytic_Data, V4 < 8.8)}



#Simulation_Data[5] = Simulation_Data[5] - log(Simulation_Data[2])
#Semianalytic_Data[5] = Semianalytic_Data[5] - log(Semianalytic_Data[2])
#Simulation_Data[6] = Simulation_Data[6] - log(Simulation_Data[2])
#Semianalytic_Data[6] = Semianalytic_Data[6] - log(Semianalytic_Data[2])
#Simulation_Data[7] = Simulation_Data[7] - log(Simulation_Data[2])
#Semianalytic_Data[7] = Semianalytic_Data[7] - log(Semianalytic_Data[2])
colnames(Simulation_Data)[4] <- "log_of_distance"
colnames(Semianalytic_Data)[4] <- "log_of_distance"
colnames(Simulation_Data)[5] <- "log_of_mean_homozygosity"
colnames(Semianalytic_Data)[5] <- "log_of_mean_homozygosity"

Semianalytic_Data_rho_1  <-subset(Semianalytic_Data, V2 ==1)
Semianalytic_Data_rho_p1  <-subset(Semianalytic_Data, V2 ==.1)
Semianalytic_Data_rho_p01  <-subset(Semianalytic_Data, V2 ==.01)

p <- ggplot() + geom_smooth(data=Semianalytic_Data_rho_1, aes(log_of_distance, log_of_mean_homozygosity), color = "green", se = FALSE)   + geom_smooth(data=Semianalytic_Data_rho_p1, aes(log_of_distance, log_of_mean_homozygosity), color = "yellow", se = FALSE) + geom_smooth(data=Semianalytic_Data_rho_p01, aes(log_of_distance, log_of_mean_homozygosity), color = "orange", se = FALSE) + 
        geom_pointrange(data=Simulation_Data, aes(x = log_of_distance, y =log_of_mean_homozygosity, ymin =V6, ymax =V7), color="black", pch = 0)  + ggtitle(paste("Alpha", alpha, "  Mu", mu)) + labs(x= "Log of Distance"
, y = "Log of Mean Homozygosity" )  


alpha_dummy <- (alpha -1)*100
file_name <- paste("log_plot_MH_varying_rho_alpha_1p", alpha_dummy, ".pdf", sep = "")
#pdf("plots.pdf")
pdf(file_name)
print(p)


dev.off()

}



