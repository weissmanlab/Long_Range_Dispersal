library(ggplot2)
library(gridExtra)
library(Hmisc)
library(dplyr)

tail_cutoff <- 5
par(mfrow=c(2,2))
x <- list()
y <- list()
p <- list(4)
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

#if(i==0)
#{mu <- .1}

#if(i==1)
#{mu <- .01}

#if(i==2)
#{mu <- .001}



Simulation_Data <-   read.table("Total_outfile_ALL.txt")
Semianalytic_Data <- read.table("formatted_semianalytic_data_ALL.txt")

Simulation_Data  <- subset(Simulation_Data, V1 == alpha)
#Simulation_Data  <-subset(Simulation_Data, V2 ==rho_inverse)
#Simulation_Data <- subset(Simulation_Data, V3 ==mu)

Semianalytic_Data <- read.table("formatted_semianalytic_data_ALL.txt")
Semianalytic_Data  <- subset(Semianalytic_Data, V1 == alpha)
#Semianalytic_Data  <-subset(Semianalytic_Data, V2 ==rho_inverse)
#Semianalytic_Data  <- subset(Semianalytic_Data, V3 ==mu)






Simulation_Data <- Simulation_Data[!(Simulation_Data$V3 == 1 ),]
Semianalytic_Data  <-Semianalytic_Data [!(Semianalytic_Data$V3 == 1 ),]

if (alpha > 1.45)
{#Simulation_Data <- Simulation_Data[!(Simulation_Data$V7 - Simulation_Data$V6  > 5),]
#Semianalytic_Data <- Semianalytic_Data[!(Simulation_Data$V7 - Simulation_Data$V6  > 5),]

Simulation_Data <- Simulation_Data[!(Simulation_Data$V4 > 10),]
Semianalytic_Data <- Semianalytic_Data[!(Semianalytic_Data$V4 > 10),]
# last distance point has extremely large error bars.  Remove this simulation data w/ insufficient trials.
}
#Simulation_Data <- Simulation_Data[!(Simulation_Data$V5   < -30),]



## convert from log probability to log odds 
Simulation_Data[5] <- log(exp(Simulation_Data[5])/(1 - exp(Simulation_Data[5]))) -log(Simulation_Data[2])
Semianalytic_Data[5] <- log(exp(Semianalytic_Data[5])/(1 - exp(Semianalytic_Data[5]))) -log(Semianalytic_Data[2])
#print(Simulation_Data[5])
Simulation_Data[6] <- log(exp(Simulation_Data[6])/(1 - exp(Simulation_Data[6]))) -log(Simulation_Data[2])
Semianalytic_Data[6] <- log(exp(Semianalytic_Data[6])/(1 - exp(Semianalytic_Data[6]))) -log(Semianalytic_Data[2])

Simulation_Data[7] <- log(exp(Simulation_Data[7])/(1 - exp(Simulation_Data[7]))) -log(Simulation_Data[2])
Semianalytic_Data[7] <- log(exp(Semianalytic_Data[7])/(1 - exp(Semianalytic_Data[7]))) -log(Semianalytic_Data[2])


#rescale log odds to collapse plots w/ different mu
Simulation_Data[4] <- Simulation_Data[4] + (1/Simulation_Data[1])*log(Simulation_Data[3])
Semianalytic_Data[4] <- Semianalytic_Data[4] + (1/Semianalytic_Data[1])*log(Semianalytic_Data[3])

Simulation_Data[5] <- Simulation_Data[5] - (1/Simulation_Data[1] - 1)*log(Simulation_Data[3])
Semianalytic_Data[5] <- Semianalytic_Data[5] - (1/Semianalytic_Data[1] - 1)*log(Semianalytic_Data[3])

Simulation_Data[6] <- Simulation_Data[6] - (1/Simulation_Data[1] - 1)*log(Simulation_Data[3])
Semianalytic_Data[6] <- Semianalytic_Data[6] - (1/Semianalytic_Data[1] - 1)*log(Semianalytic_Data[3])

Simulation_Data[7] <- Simulation_Data[7] - (1/Simulation_Data[1] - 1)*log(Simulation_Data[3])
Semianalytic_Data[7] <- Semianalytic_Data[7] - (1/Semianalytic_Data[1] - 1)*log(Semianalytic_Data[3])
#replace mean homozygosity/probability of identity with log of odds of identiy + log(rho)


#take tail for power law estimation
Tail_Sum <- 0
Simulation_Data_Tail <- subset(Simulation_Data, V4 >= tail_cutoff)
Tail_Sum <- sum(exp(Simulation_Data_Tail$V5)) # interpret probabilites as frequencies.  sum of frequencies is analogous to number of trials/events in the tail of a histogram divided by the total number of trials
Expected_log_dist_from_xmin <- sum(exp(Simulation_Data_Tail$V5)*(Simulation_Data_Tail$V4 -log(exp(tail_cutoff) )))/Tail_Sum

tail_estimate <- 1 + 1/Expected_log_dist_from_xmin
alpha_estimate <- tail_estimate -1
print(tail_estimate) 
#See "Power Law Distributions in Empirical Data" by Clauset, Shalizi and Newman for an explanation of this maximum likelihood algorithm for power law tails.  The formula they give in section three to estimate tails from discrete data in a histogram is alpha = 1+ n[Sum(ln(xi/xmin))]^-1, where the sum is taken over all trials found in a bin greater than or equal to the tail cutoff xmin.  Since we have frequencies rather that a certain number of trials at each distance, we need to re-express this formula in terms of frequencies.  alpha = 1+ n[Sum(n_i ln(xi/xmin))]^-1 where the sum is now taken over all x.  This is equivlaent to 1+ f_tail[Sum(f_i ln(xi/xmin))]^-1, where f_i is the fraction of points found at a particular distance and f_tail = sum(f_i) is the fraction of points found in the tail.  f_i/f_tail = n_i/n, and the estimators are equvialent. 


#Simulation_Data[6] = log(exp(Simulation_Data[6])(1 - exp(Simulation_Data[6]))) -#log(Simulation_Data[2])
#Simulation_Data[7] = log(exp(Simulation_Data[7])(1 - exp(Simulation_Data[7]))) -log(Simulation_Data[2])
#Semianalytic_Data[5] = Semianalytic_Data[5] - log(Semianalytic_Data[2])
#Simulation_Data[6] = Simulation_Data[6] - log(Simulation_Data[2])
#Semianalytic_Data[5] = Semianalytic_Data[5] - log(Semianalytic_Data[2])
#Simulation_Data[6] = Simulation_Data[6] - log(Simulation_Data[2])
#Semianalytic_Data[7] = Semianalytic_Data[7] - log(Semianalytic_Data[2])
colnames(Simulation_Data)[4] <- "log_of_distance_rescaled"
colnames(Semianalytic_Data)[4] <- "log_of_distance_rescaled"
colnames(Simulation_Data)[5] <- "log_odds_of_identity_rescaled"
colnames(Semianalytic_Data)[5] <- "log_odds_of_identity_rescaled"





Rho <- as.factor(1/Simulation_Data[, 2])
Mu <- as.factor(Simulation_Data[, 3])

#colnames(Simulation_Data["V2"]) = "Rho_Inverse"
#colnames(Simulation_Data["V3"]) = "Mu"
#Simulation_Data["COLOR"] <- as.factor(-(log10(Simulation_Data[,3]) -14))
#Simulation_Data["SHAPE"] <- as.factor(Simulation_Data[,2]) #as.factor(-log10(Simulation_Data[,2]) +1)


p[[j]] <- ggplot(data=Simulation_Data, aes(x = log_of_distance_rescaled, y =log_odds_of_identity_rescaled, ymin =V6, ymax =V7, colour = Mu, shape = Rho)) + geom_point() +geom_pointrange() + ggtitle(paste("alpha", alpha)) + geom_smooth(data=Semianalytic_Data, colour = "brown") 


#p[[k]]

plot(p[[j]])


#k <- k + 1

#basic <- ggplot(mtcars, aes(wt, mpg, colour = factor(cyl), shape = factor(vs) )) + geom_point()



#}


}

#multiplot(p)
#do.call(grid.arrange,p)

####################


