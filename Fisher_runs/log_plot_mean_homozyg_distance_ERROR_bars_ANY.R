library(ggplot2)
library(gridExtra)
library(Hmisc)
library(dplyr)
alpha <- 2.05
tail_cutoff <- 111 # on natural log scale
#par(mfrow=c(3,3))
x <- list()
y <- list()
p <- list()
k <- 1
for(i in 0:2){for(j in 1:3){

rho_inverse <- 10^(-i)
mu <- 10^(-j)
#dummy_slope_Data <-   read.table("dummy_slope.txt")
Simulation_Data <-   read.table("ALL_Fisher_runs_MH.txt")
Simulation_Data  <- subset(Simulation_Data, V1 == alpha)
Simulation_Data  <-subset(Simulation_Data, V2 ==rho_inverse)
Simulation_Data <- subset(Simulation_Data, V3 ==mu)


Simulation_Data[,4] <- log(Simulation_Data[,4])
Simulation_Data[,5] <- log(Simulation_Data[,5])
Simulation_Data[,6] <- log(Simulation_Data[,6])
Simulation_Data[,7] <- log(Simulation_Data[,7])

#Semianalytic_Data <- read.table("formatted_semianalytic_data_ALL.txt")
#Semianalytic_Data  <- subset(Semianalytic_Data, V1 == alpha)
#Semianalytic_Data  <-subset(Semianalytic_Data, V2 ==rho_inverse)
#Semianalytic_Data  <- subset(Semianalytic_Data, V3 ==mu)

#All_Data <- rbind(Simulation_Data,Semianalytic_Data)
#All_Data$labels[(length(Simulation_Data[,1])+1):length(All_Data[,1])] <- 'ni2'
#All_Data$labels[1:length(Simulation_Data[,1])] <- 'ni1'

#All_Data <- All_Data[!(All_Data$V7 - All_Data$V6  > 5),]
Simulation_Data <- Simulation_Data[!(Simulation_Data$V7 - Simulation_Data$V6  > 5),]
#Semianalytic_Data <- Semianalytic_Data[!(Semianalytic_Data$V7 - Semianalytic_Data$V6  > 5),]


#take tail for power law estimation
#Tail_Sum <- 0
#Dummy_Data_Tail <- subset(dummy_slope_Data, V1 >= tail_cutoff)
#Tail_Sum <- sum(exp(Dummy_Data_Tail$V2)) # interpret probabilites as frequencies.  sum of frequencies is analogous to number of trials/events in the tail of a histogram divided by the total number of trials
#Expected_log_dist_from_xmin <- sum(exp(Dummy_Data_Tail$V2)*(Dummy_Data_Tail$V1 -log(exp(tail_cutoff) )))/Tail_Sum

#tail_estimate <- 1 + 1/Expected_log_dist_from_xmin
#alpha_estimate <- tail_estimate -1
#print(tail_estimate) 
#See "Power Law Distributions in Empirical Data" by Clauset, Shalizi and Newman for an explanation of this maximum likelihood algorithm for power law tails.  The formula they give in section three to estimate tails from discrete data in a histogram is alpha = 1+ n[Sum(ln(xi/xmin))]^-1, where the sum is taken over all trials found in a bin greater than or equal to the tail cutoff xmin.  Since we have frequencies rather that a certain number of trials at each distance, we need to re-express this formula in terms of frequencies.  alpha = 1+ n[Sum(n_i ln(xi/xmin))]^-1 where the sum is now taken over all x.  This is equivlaent to 1+ f_tail[Sum(f_i ln(xi/xmin))]^-1, where f_i is the fraction of points found at a particular distance and f_tail = sum(f_i) is the fraction of points found in the tail.  f_i/f_tail = n_i/n, and the estimators are equvialent. 





colnames(Simulation_Data)[4] <- "log_of_distance"
#colnames(Semianalytic_Data)[4] <- "log_of_distance"
colnames(Simulation_Data)[5] <- "log_of_mean_homozygosity"
#colnames(Semianalytic_Data)[5] <- "log_of_mean_homozygosity"


p[[k]] <- ggplot() + geom_pointrange(data=Simulation_Data, aes(x = log_of_distance, y =log_of_mean_homozygosity, ymin =V6, ymax =V7), color="blue", pch = 0)  + ggtitle(paste("alpha", alpha, "rho", 1/rho_inverse, "mu", mu))
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
