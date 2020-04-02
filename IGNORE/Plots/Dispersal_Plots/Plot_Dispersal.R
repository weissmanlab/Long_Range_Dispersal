library(ggplot2)
library(libstableR)
library(argparse)
library(mdatools)

scale_parameter <- 250
#sigma <- c/4
alpha1 <- 2
alpha2 <- 1.45
time1 <- .01
time2 <- .1
time3 <- 1
time4 <- 10
time5 <- 100
initial_distance <- 1000

#Make dataframe for Gaussinan profile at various times
#Make dataframe for Levy alpha stable profile at various times
#plot and compare to get a sense for long range vs short range dispersal

Gaussian_Dispersal_Data <- data.frame(matrix(0, 200,  6))
Gaussian_Dispersal_Data[,1] <- alpha1

Stable_Dispersal_Data <- data.frame(matrix(0, 200,  6))
Stable_Dispersal_Data[,1] <- alpha2
for(i in 1:200)
{ 
  index <- i + 1
  log_dist <- i/10
  dist <- exp(log_dist)
  Gaussian_Dispersal_Data[index, 2] <- log_dist
  Stable_Dispersal_Data[index, 2] <- log_dist 
  time <- 1
  pars1 <- c(alpha1, 0, scale_parameter*((time)^(1/alpha1)), 0)
  pars2 <- c(alpha2, 0, scale_parameter*((time)^(1/alpha2)), 0)
  Stable_Dispersal_Data[index, 3] <- log(2*stable_pdf(dist , pars2))
  Gaussian_Dispersal_Data[index, 3] <- log(2*dnorm(dist , 0, .5*scale_parameter*(time)^(1/alpha1)))
  
  time <- 10
  pars1 <- c(alpha1, 0, scale_parameter*((time)^(1/alpha1)), 0)
  pars2 <- c(alpha2, 0, scale_parameter*((time)^(1/alpha2)), 0)
  Stable_Dispersal_Data[index, 4] <- log(2*stable_pdf(dist , pars2))
   Gaussian_Dispersal_Data[index, 4] <- log(2*dnorm(dist , 0, .5*scale_parameter*(time)^(1/alpha1)))
  
  time <- 100
  pars1 <- c(alpha1, 0, scale_parameter*((time)^(1/alpha1)), 0)
  pars2 <- c(alpha2, 0, scale_parameter*((time)^(1/alpha2)), 0)
  Stable_Dispersal_Data[index, 5] <- log(2*stable_pdf(dist , pars2))
   Gaussian_Dispersal_Data[index, 5] <- log(2*dnorm(dist , 0, .5*scale_parameter*(time)^(1/alpha1)))
  
  time <- 1000
  pars1 <- c(alpha1, 0, scale_parameter*((time)^(1/alpha1)), 0)
  pars2 <- c(alpha2, 0, scale_parameter*((time)^(1/alpha2)), 0)
  Stable_Dispersal_Data[index, 6] <- log(2*stable_pdf(dist , pars2))
   Gaussian_Dispersal_Data[index, 6] <- log(2*dnorm(dist , 0, .5*scale_parameter*(time)^(1/alpha1)))
	

}

#print(Gaussian_Dispersal_Data)
#p <- ggplot() + geom_smooth(data=Gaussian_Dispersal_Data, aes(x = Gaussian_Dispersal_Data[, 2], y = Gaussian_Dispersal_Data[, 4], color = "Gaussian"), se = FALSE) + geom_smooth(data=Stable_Dispersal_Data, aes(x = Stable_Dispersal_Data[, 2], y = Stable_Dispersal_Data[, 3], color = "Stable"), se = FALSE)
pdf("Levy_alpha_stable_dispersal_varying_time.pdf")

p <- ggplot() + geom_smooth(data=Stable_Dispersal_Data, aes(x = Stable_Dispersal_Data[, 2], y = Stable_Dispersal_Data[, 3], color = "time = 1"), se = FALSE) + geom_smooth(data=Stable_Dispersal_Data, aes(x = Stable_Dispersal_Data[, 2], y = Stable_Dispersal_Data[, 4], color = "time = 10"), se = FALSE) + geom_smooth(data=Stable_Dispersal_Data, aes(x = Stable_Dispersal_Data[, 2], y = Stable_Dispersal_Data[, 5], color = "time = 100"), se = FALSE) + geom_smooth(data=Stable_Dispersal_Data, aes(x = Stable_Dispersal_Data[, 2], y = Stable_Dispersal_Data[, 6], color = "time = 1000"), se = FALSE)  + scale_colour_manual("", breaks = c("time = 1", "time = 10", "time = 100", "time = 1000"), values = c("red", "orange", "brown", "black"))+labs( x = "Log of Distance", y= "Log of Probability Density") +ggtitle("Levy Flight Dispersal")


print(p)

dev.off()


#Gaussian_Dispersal_Data <- subset( Gaussian_Dispersal_Data, Gaussian_Dispersal_Data[,2] < 10 )

pdf("Gaussian_dispersal_varying_time.pdf")

p <- ggplot() + geom_smooth(data=Gaussian_Dispersal_Data, aes(x = Gaussian_Dispersal_Data[, 2], y = Gaussian_Dispersal_Data[, 3], color = "time = 1"), se = FALSE) + geom_smooth(data=Gaussian_Dispersal_Data, aes(x = Gaussian_Dispersal_Data[, 2], y = Gaussian_Dispersal_Data[, 4], color = "time = 10"), se = FALSE) + geom_smooth(data=Gaussian_Dispersal_Data, aes(x = Gaussian_Dispersal_Data[, 2], y = Gaussian_Dispersal_Data[, 5], color = "time = 100"), se = FALSE) + geom_smooth(data=Gaussian_Dispersal_Data, aes(x = Gaussian_Dispersal_Data[, 2], y = Gaussian_Dispersal_Data[, 6], color = "time = 1000"), se = FALSE)  + scale_colour_manual("", breaks = c("time = 1", "time = 10", "time = 100", "time = 1000"), values = c("red", "orange", "brown", "black"))+labs( x = "Log of Distance", y= "Log of Probability Density") +ggtitle("Gaussian Dispersal")


print(p)

dev.off()





#print(dummy)