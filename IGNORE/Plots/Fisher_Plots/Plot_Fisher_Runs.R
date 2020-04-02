library(ggplot2)
library(gridExtra)
library(Hmisc)
library(dplyr)
library(mosaic)
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
c <- 5 # 250
sigma <- c
#dummy_slope_Data <-   read.table("dummy_slope.txt")
Simulation_Data <-   read.table("ALL_Fisher_runs_MH.txt")
Simulation_Data  <- subset(Simulation_Data, V1 == alpha)
Simulation_Data  <-subset(Simulation_Data, V2 ==rho_inverse)
Simulation_Data <- subset(Simulation_Data, V3 ==mu)
if(mu == .001) {Simulation_Data <- subset(Simulation_Data, V4 < exp(9) + 1)
	xlimit <- exp(9) + 1
	}
if(mu == .01) {Simulation_Data <- subset(Simulation_Data, V4 <exp(8) +1)
	xlimit <- exp(8) +1 
	}
if(mu == .1) {Simulation_Data <- subset(Simulation_Data, V4 < exp(7) + 1)
	xlimit <- exp(7) + 1
	}

#print(mu)
Gaussian_dispersal_prediction <- (1/(1 + (2*sigma*sqrt(2*mu)/rho_inverse )))*exp(-(sqrt(2*mu)/sigma) *Simulation_Data[,4])


#Gaussian_dispersal_prediction <- exp(-(sqrt(2*mu)/sigma) *Simulation_Data[,4])
#Gaussian_dispersal_prediction <- exp(- .01*Simulation_Data[,4])
#print(Gaussian_dispersal_prediction)

analytic_data <- data.frame(Simulation_Data[, 1], Simulation_Data[, 2], Simulation_Data[, 3], Simulation_Data[, 4], Gaussian_dispersal_prediction, Gaussian_dispersal_prediction, Gaussian_dispersal_prediction)
#print(analytic_data)
#print(Simulation_Data[,4])

rescaled_Simulation_Data <- Simulation_Data

exponential_fit <-  fitModel(V5 ~ (a^2)*exp(-(b^2)*V4), data =Simulation_Data)


exponential_function <- data.frame(matrix(0, ncol = 2, nrow = xlimit))
#exponential_function[1, 1] <- 0
for( count in 1:xlimit){exponential_function[count, 1] = count-1}



exponential_function[,2] <- exponential_fit(exponential_function[,1])

#print(exponential_function)


#print(exponential_fit(5))

p[[k]] <- ggplot() + geom_pointrange(data=Simulation_Data, aes(x = V4, y =V5, ymin =V6, ymax =V7, color="Simulation Values"), pch = 0)  + ggtitle(paste("Alpha", alpha, "Rho", 1/rho_inverse, "Mu", mu)) + labs(x= "Distance", y =  "Mean Homozygosity" ) + geom_smooth(data=exponential_function, aes(x = exponential_function[,1], y =exponential_function[,2], color="Exponential Fit", se = FALSE), pch = 0)  + scale_colour_manual("", breaks = c("Simulation Values", "Exponential Fit"), values = c("red", "blue")) 

#p[[k]] <-  ggplot() + stat_function(fun =x exponential_fit) + xlim(0,3000)

#+ geom_pointrange(data = analytic_data, aes(x = analytic_data[,4], y =analytic_data[,5], ymin =analytic_data[,6], ymax =analytic_data[,7]), color = "brown")

#p[[k]] <- ggplot() + geom_pointrange(data = analytic_data, aes(x = analytic_data[,4], y =analytic_data[,5], ymin =analytic_data[,6], ymax =analytic_data[,7]), color = "brown")


 #p <- p + geom_point()
#print(p)
#p <- ggplot(mtcars, aes(wt, mpg))
#p <- p + geom_point()
if(mu == .1){ mu_dummy <- "p1"}
if(mu == .01){ mu_dummy <- "p01"}
if(mu == .001){ mu_dummy <- "p001"}

pdf(paste("Fisher_plot_rho", 1/rho_inverse, "mu_", mu_dummy, ".pdf", sep = ""))

print(p[[k]])



dev.off()
k <- k + 1

###############
#dist <- ts(dummy)
#plot.ts( dist)
#plot(dummy)
#abline(h = 4, v = 4, col = "gray60")
#abline(v = 4, col = "gray60")
#abline(a = 9.8, b = -2.95)
 
 
 


}}

#do.call(grid.arrange,p)



####################
