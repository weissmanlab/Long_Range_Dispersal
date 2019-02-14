library(Hmisc)
library(dplyr)
alpha <- 1.85
par(mfrow=c(3,3))

for(i in 0:2){for(j in 1:3){

rho_inverse <- 10^(-i)
mu <- 10^(-j)

Simulation_Data <-   read.table("Total_outfile_ALL.txt")
Simulation_Data  <- subset(Simulation_Data, V1 == alpha)
Simulation_Data  <-subset(Simulation_Data, V2 ==rho_inverse)
Simulation_Data <- subset(Simulation_Data, V3 ==mu)

Semianalytic_Data <- read.table("formatted_semianalytic_data_ALL.txt")
Semianalytic_Data  <- subset(Semianalytic_Data, V1 == alpha)
Semianalytic_Data  <-subset(Semianalytic_Data, V2 ==rho_inverse)
Semianalytic_Data <- subset(Semianalytic_Data, V3 ==mu)



#dist <- ts(dummy)
#plot.ts( dist)
#plot(dummy)
#abline(h = 4, v = 4, col = "gray60")
#abline(v = 4, col = "gray60")
#abline(a = 9.8, b = -2.95)
 
 
 
 distance <- as.vector(Semianalytic_Data[,4])
 class(distance)
 avg <- as.vector(Semianalytic_Data[,5])
 class(avg)
 upper <- Semianalytic_Data[,6]
 lower <- Semianalytic_Data[,7]
 plot(distance, avg, col = "red", ylim=range(c(lower, upper)),
    pch=19, xlab = "Log of Distance", ylab = "Log of Mean Homozygosity", main = paste("alpha", alpha, "rho", 1/rho_inverse, "mu", mu))

#arrows(distance, lower, distance, upper, length=0.05, angle=90, code=3)

 
  
distance <- as.vector(Simulation_Data[,4])
class(distance)
avg <- as.vector(Simulation_Data[,5])
class(avg)
upper <- Simulation_Data[,6]
lower <- Simulation_Data[,7]
points(distance, avg, col ="blue" ,ylim=range(c(lower, upper)), pch=19)
arrows(distance, lower, distance, upper, length=0.05, angle=90, code=3, col = "blue")



legend( x= "bottomleft", y=0.99, 
        legend=c("Simulations","Numerics"), 
        col=c("blue",  "red"),   
        pch=c(19, 19))

 
}}
