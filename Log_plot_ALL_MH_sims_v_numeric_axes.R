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


Simulation_Data <-   read.table("ALL_runs_MH.txt")
Semianalytic_Data <-read.table("numeric_MH_LOG_SCALE_ALL_but_zeros.txt")

Simulation_Data <- Simulation_Data[!(Simulation_Data$V3 == 1 ),]

Simulation_Data <- subset(Simulation_Data, V1 < 2) # alpha > 2 has no numerics
#convert simulation data to log scale
Simulation_Data <- subset(Simulation_Data, V4 > 0)
# can't put zeros on logs scale, they diverge
Simulation_Data[,4] <-log(Simulation_Data[,4])
Simulation_Data[,5] <-log(Simulation_Data[,5])
Simulation_Data[,6] <-log(Simulation_Data[,6])
Simulation_Data[,7] <-log(Simulation_Data[,7])

Simulation_Data <- subset(Simulation_Data, V3 < .9)
Semianalytic_Data  <- subset(Semianalytic_Data, V3 <.9)

Simulation_Data <- subset(Simulation_Data, V4 < 10.9)
Semianalytic_Data  <- subset(Semianalytic_Data, V4 <10.9)


## now we must sort both dataframes so that all the entries are in the same order
Simulation_Data <- Simulation_Data[order(Simulation_Data[, 4], Simulation_Data[, 3], Simulation_Data[, 2], Simulation_Data[, 1]), ]


Semianalytic_Data <- Semianalytic_Data[order(Semianalytic_Data[,4], Semianalytic_Data[,3], Semianalytic_Data[,2], Semianalytic_Data[,1] ),]

#Semianalytic_Data[order(Semianalytic_Data[,4]),]
#print(Simulation_Data - Semianalytic_Data)
#print( Semianalytic_Data)
#print(Simulation_Data)


numeric_axis <- Semianalytic_Data[, 5]
simulation_axis <- Simulation_Data[, 5]
lower_CI <- Simulation_Data[, 6]
upper_CI <- Simulation_Data[, 7]

#rbind(numeric_axis, simulation_axis)

df <- data.frame(numeric_axis, simulation_axis, lower_CI, upper_CI)



p <- ggplot() + geom_pointrange(data=df, aes(x = numeric_axis, y =simulation_axis, ymin =lower_CI, ymax =upper_CI), color="blue", pch = 0)  + ggtitle(paste("Log of Mean Homozygosity Simulations vs Numerics"))

pdf("Log_plot_MH_sims_v_numerics.pdf")
print(p)
dev.off()