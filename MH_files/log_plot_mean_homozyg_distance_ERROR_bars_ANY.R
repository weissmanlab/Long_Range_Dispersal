library(ggplot2)
library(gridExtra)
library(Hmisc)
library(dplyr)
alpha <- 1.25
#par(mfrow=c(3,3))
x <- list()
y <- list()
p <- list()
k <- 1
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
Semianalytic_Data  <- subset(Semianalytic_Data, V3 ==mu)

All_Data <- rbind(Simulation_Data,Semianalytic_Data)
All_Data$labels[(length(Simulation_Data[,1])+1):length(All_Data[,1])] <- 'ni2'
All_Data$labels[1:length(Simulation_Data[,1])] <- 'ni1'

All_Data <- All_Data[!(All_Data$V7 - All_Data$V6  > 5),]
Simulation_Data <- Simulation_Data[!(Simulation_Data$V7 - Simulation_Data$V6  > 5),]
Semianalytic_Data <- Semianalytic_Data[!(Semianalytic_Data$V7 - Semianalytic_Data$V6  > 5),]


colnames(Simulation_Data)[4] <- "log_of_distance"
colnames(Semianalytic_Data)[4] <- "log_of_distance"
colnames(Simulation_Data)[5] <- "log_of_mean_homozygosity"
colnames(Semianalytic_Data)[5] <- "log_of_mean_homozygosity"


p[[k]] <- ggplot() + geom_point(data=Semianalytic_Data, aes(log_of_distance, log_of_mean_homozygosity), color = "red") + 
        geom_pointrange(data=Simulation_Data, aes(x = log_of_distance, y =log_of_mean_homozygosity, ymin =V6, ymax =V7), color="blue", pch = 0)  + ggtitle(paste("alpha", alpha, "rho", 1/rho_inverse, "mu", mu))
k <- k + 1
#p <- p + geom_point()
#print(p)
#p <- ggplot(mtcars, aes(wt, mpg))
#p <- p + geom_point()
#print(p)


###############
#dist <- ts(dummy)
#plot.ts( dist)
#plot(dummy)
#abline(h = 4, v = 4, col = "gray60")
#abline(v = 4, col = "gray60")
#abline(a = 9.8, b = -2.95)
 
 
 if(FALSE) {

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
#par(new=TRUE)
points(distance, avg, col ="blue" ,ylim=range(c(lower, upper)), pch=19)
arrows(distance, lower, distance, upper, length=0.05, angle=90, code=3, col = "blue")



legend( x= "bottomleft", y=0.99, 
        legend=c("Simulations","Numerics"), 
        col=c("blue",  "red"),   
        pch=c(19, 19))
}
 



}}

do.call(grid.arrange,p)



####################