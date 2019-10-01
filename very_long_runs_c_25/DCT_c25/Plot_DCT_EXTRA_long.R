library(libstableR)
library(argparse)
library(mdatools)
library(ggplot2)
library(gridExtra)
library(Hmisc)
library(dplyr)
library(ggpmisc)
Final_Time <-100000
Coalescence_Data_dummy <-   read.table("ALL_runs_DCT.txt")
Coalescence_Data_dummy <- head(Coalescence_Data_dummy, Final_Time)
for(j in 1:9){for(ZZ in 0:0){
if(j==1)
{alpha <- 1.25}
if(j==2)
{alpha <- 1.45}
if(j==3)
{alpha <- 1.65}
if(j==4)
{alpha <- 1.85}
if(j==5)
{alpha <- 2.05}
if(j==6)
{alpha <- .25}
if(j==7)
{alpha <- .5}
if(j==8)
{alpha <- .75}
if(j==9)
{alpha <- 1}


#alpha <- 1.65
distance <- 0 #exp(ZZ)
distance_lower_bound <- .9*distance - .001
distance_upper_bound <- 1.1*distance + .001




Coalescence_Data_plot_ALL <-  data.frame(matrix(0, Final_Time, 2))
Numeric_approx_Data_plot_ALL <-  data.frame(matrix(0, Final_Time, 2))
Asymptotic_approx_Data_plot_ALL <-  data.frame(matrix(0, Final_Time, 2))
for( i in 1:Final_Time )
{   Numeric_approx_Data_plot_ALL[i,1] <- log(i)
	Asymptotic_approx_Data_plot_ALL[i,1] <- log(i)
	time_dummy <- i 
	 
	 scale_parameter <- 25.0
	 generalized_D <- (scale_parameter^alpha)/2
	 pars <- c(alpha, 0, scale_parameter*((time_dummy)^(1/alpha)), 0) 
     #log_stable_dist <- log(stable_pdf(distance, pars))	
	#Numeric_approx_Data_plot_ALL[i,2] <-  log_stable_dist + log(1)
	if(alpha < 1){Asymptotic_approx_Data_plot_ALL[i,2] <-  log((gamma(1+1/alpha)*(2*generalized_D*time_dummy)^(-1/alpha))/(pi))}
	if(alpha == 1){Asymptotic_approx_Data_plot_ALL[i,2] <-  log(2*pi*generalized_D /(time_dummy* log(scale_parameter*time_dummy)^2))}
	if(alpha > 1){Asymptotic_approx_Data_plot_ALL[i,2] <-  log(scale_parameter*(alpha -1)*sin(pi/alpha)*(time_dummy)^(1/alpha -2))}

	
	#Numeric_approx_Data_plot_ALL[i,3] <- log_stable_dist + log(1/10)
	#Numeric_approx_Data_plot_ALL[i,4] <- log_stable_dist + log(1/100)
	
	}



for( Q in 1:1)
{
rho_inverse <- 1 #10^(-1*(Q-1))
rho <- 1/rho_inverse
#print(rho_inverse)

Coalescence_Data <- Coalescence_Data_dummy



Coalescence_Data <- subset(Coalescence_Data, V1 == alpha)
Coalescence_Data <- subset(Coalescence_Data, V2 == 1/rho)
Coalescence_Data <- subset(Coalescence_Data, V3 > distance_lower_bound)
Coalescence_Data <- subset(Coalescence_Data, V3 < distance_upper_bound)

Coalescence_Data_Transposed <- t(Coalescence_Data)
#Coalescence_Data_Transposed <- Coalescence_Data_Transposed[-c(1:3), ]


Coalescence_Data_plot <-  data.frame(matrix(0, Final_Time, 2))
for( i in 1:Final_Time )
{Coalescence_Data_plot[i, 1] <- log(i)
 Coalescence_Data_plot[i, 2] <- log(Coalescence_Data_Transposed[i + 3, 1])	

 	
}

#print(Coalescence_Data_plot[, 2])
#normalization_check <- sum(exp(Coalescence_Data_plot[, 2]))
#normalization_check <- sum(Coalescence_Data_plot[, 2])
#print(normalization_check)
Coalescence_Data_plot_ALL[, 1] <-Coalescence_Data_plot[, 1]
Coalescence_Data_plot_ALL[, Q+ 1] <-Coalescence_Data_plot[, 2]
}


#print(Coalescence_Data_plot)

#Coalescence_Data_plot <- subset(Coalescence_Data_plot

alpha_dummy <- 100*(alpha -1)
if(j < 5){
#print(paste("Log_Plot_DCT_LONG_alpha_1p", alpha_dummy, "distance_e", trunc(log(distance)), ".pdf",sep=""))
#pdf(paste("Log_Plot_DCT_LONG_alpha_1p", alpha_dummy, "distance_e", trunc(log(distance)), ".pdf",sep=""))
print(paste("Log_Plot_DCT_LONG_alpha_1p", alpha_dummy, "distance_0" ,".pdf",sep=""))
pdf(paste("Log_Plot_DCT_LONG_alpha_1p", alpha_dummy, "distance_0", ".pdf",sep=""))

}
if(j == 5){
#print(paste("Log_Plot_DCT_LONG_alpha_2p05", "distance_e", trunc(log(distance)), ".pdf",sep=""))
#pdf(paste("Log_Plot_DCT_LONG_alpha_2p05", "distance_e", trunc(log(distance)), ".pdf",sep=""))
print(paste("Log_Plot_DCT_LONG_alpha_2p05", "distance_0", ".pdf",sep=""))
pdf(paste("Log_Plot_DCT_LONG_alpha_2p05", "distance_0", ".pdf",sep=""))
}
if(j == 6){
#print(paste("Log_Plot_DCT_LONG_alpha_0p25", "distance_e", trunc(log(distance)), ".pdf",sep=""))
#pdf(paste("Log_Plot_DCT_LONG_alpha_0p25", "distance_e", trunc(log(distance)), ".pdf",sep=""))
print(paste("Log_Plot_DCT_LONG_alpha_0p25", "distance_0", ".pdf",sep=""))
pdf(paste("Log_Plot_DCT_LONG_alpha_0p25", "distance_0", ".pdf",sep=""))
}
if(j == 7){
#print(paste("Log_Plot_DCT_LONG_alpha_0p5", "distance_e", trunc(log(distance)), ".pdf",sep=""))
#pdf(paste("Log_Plot_DCT_LONG_alpha_0p5", "distance_e", trunc(log(distance)), ".pdf",sep=""))
print(paste("Log_Plot_DCT_LONG_alpha_0p5", "distance_0", ".pdf",sep=""))
pdf(paste("Log_Plot_DCT_LONG_alpha_0p5", "distance_0", ".pdf",sep=""))

}
if(j == 8){
#print(paste("Log_Plot_DCT_LONG_alpha_0p75", "distance_e", trunc(log(distance)), ".pdf",sep=""))
#pdf(paste("Log_Plot_DCT_LONG_alpha_0p75", "distance_e", trunc(log(distance)), ".pdf",sep=""))
print(paste("Log_Plot_DCT_LONG_alpha_0p75", "distance_0", ".pdf",sep=""))
pdf(paste("Log_Plot_DCT_LONG_alpha_0p75", "distance_0", ".pdf",sep=""))}
if(j == 9){
#print(paste("Log_Plot_DCT_LONG_alpha_1", "distance_e", trunc(log(distance)), ".pdf",sep=""))
#pdf(paste("Log_Plot_DCT_LONG_alpha_1", "distance_e", trunc(log(distance)), ".pdf",sep=""))
print(paste("Log_Plot_DCT_LONG_alpha_1", "distance_0", ".pdf",sep=""))
pdf(paste("Log_Plot_DCT_LONG_alpha_1", "distance_0", ".pdf",sep=""))
}


#p <- ggplot() + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X2, color = "simulated rho = 1")) + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X3, color = "simulated rho = 10"))   + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X4, color = "simulated rho = 100"))  + geom_point() + geom_smooth(data=Numeric_approx_Data_plot_ALL, aes(x = X1, y =X2, color = "numeric rho = 1", se = FALSE), linetype="dashed") + geom_smooth(data=Numeric_approx_Data_plot_ALL, aes(x = X1, y =X3, color = "numeric rho = 10", se = FALSE), linetype="dashed")   + geom_smooth(data=Numeric_approx_Data_plot_ALL, aes(x = X1, y =X4, color = "numeric rho = 100", se = FALSE), linetype="dashed")  + geom_point() + geom_point() + labs( x = "Time", y ="Log Dist of Coalescence Times") + ggtitle(paste("Alpha", alpha, "Distance e^", log(distance)))#+ labs( x = "Time", y ="Log Dist of Coalescence Times") + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X12, color = "init dist 11")) 



lm_eqn <- function(Coalescence_Data_plot_AL){
    m <- lm(y ~ x, Coalescence_Data_plot_AL);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
         list(a = format(unname(coef(m)[1]), digits = 2),
              b = format(unname(coef(m)[2]), digits = 2),
             r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
}

#my.formula = y~x

#do.call(data.frame,lapply(Coalescence_Data_plot_ALL, function(x) replace(x, is.infinite(x),NA)))
#Coalescence_Data_plot_ALL <- Coalescence_Data_plot_ALL[!is.infinite(rowSums(Coalescence_Data_plot_ALL)),]


#Coalescence_Data_plot_ALL <- na.omit(Coalescence_Data_plot_ALL)
#fit <- lm(Coalescence_Data_plot_ALL[, 2] ~ Coalescence_Data_plot_ALL[, 1], data = Coalescence_Data_plot_ALL, na.action = na.omit) 
#slope <- fit$coef[2]
#if(alpha <=1){print(-1/slope)}
#if(alpha >1){print(1/(slope +2))}



p <- ggplot() + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X2, color = "simulated rho = 1")) +   geom_point() + geom_smooth(data=Asymptotic_approx_Data_plot_ALL, aes(x = X1, y =X2, color = "asymptotic", se = FALSE), linetype="dashed")  + geom_point() + geom_point() + labs( x = "Log Time", y ="Log Dist of Coalescence Times") + ggtitle(paste("Alpha", alpha, "Distance 0")) #+ labs( x = "Time", y ="Log Dist of Coalescence Times") + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X12, color = "init dist 11")) 

if(alpha == 1)
{ p <- ggplot() + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X2, color = "simulated rho = 1")) +   geom_point() + geom_smooth(data=Asymptotic_approx_Data_plot_ALL, aes(x = X1, y =X2, color = "Kernel, NOT asymptotic", se = FALSE), linetype="dashed")  + geom_point() + geom_point() + labs( x = "Log Time", y ="Log Dist of Coalescence Times") + ggtitle(paste("Alpha", alpha, "Distance 0")) #+ labs( x = "Time", y ="Log Dist of Coalescence Times") + geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X12, color = "init dist 11")) 

}

#geom_smooth(data=Coalescence_Data_plot_ALL, aes(x = X1, y =X2, color = "simulated rho = 1"), method = "lm", se=FALSE, color="black", formula = my.formula,linetype="dashed") +




print(p)

dev.off()
#print(Coalescence_Data_plot_ALL[,4])
}}
