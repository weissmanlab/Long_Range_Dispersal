library(libstableR)
library(argparse)
library(mdatools)
library(ggplot2)
library(gridExtra)
library(Hmisc)
library(dplyr)

starting_separation <- 5
  mu <- .1
  xbar <- 1



pdf(paste("Log_Plot_DCT_varying_alpha_large_distance_limit", "_mu_", mu,"_xbar_val_", xbar, "starting_separation", starting_separation, "_xbar", ".pdf",sep=""))
p <- ggplot()
for(j in 1:3){
if(j==1)
{alpha <- .5}
if(j==2)
{alpha <- 1.5}
if(j==3)
{alpha <- 2}

 distance <- starting_separation*xbar
scale_parameter <- xbar*mu^(1/alpha)

#alpha <- 1.65

distance_lower_bound <- .9*distance - .001
distance_upper_bound <- 1.1*distance + .001


Final_Time <-1000


Numeric_approx_Data_plot_ALL <-  data.frame(matrix(0, Final_Time, 3))
for( i in 1:Final_Time )
{   Numeric_approx_Data_plot_ALL[i,1] <- i
	
	time_dummy <- i 
	 
	 	 ybar <- scale_parameter*((time_dummy)^(1/alpha)) 
	 
	 pars <- c(alpha, 0, ybar, 0) 
     log_stable_dist <- stable_pdf(distance, pars)
	Numeric_approx_Data_plot_ALL[i,2] <-  log_stable_dist 
	Numeric_approx_Data_plot_ALL[i,3] <-  alpha
	
		}
Numeric_approx_Data_plot_ALL$X3 <- factor(Numeric_approx_Data_plot_ALL$X3)


colnames(Numeric_approx_Data_plot_ALL)[3] <- "Alpha" 



p <- p + geom_smooth(data=Numeric_approx_Data_plot_ALL, aes(x = X1, y =X2, color = Alpha, se = FALSE), linetype="dashed") + geom_point() + labs( x = "Time", y ="Dist of Coalescence Times p(t)") + ggtitle(paste("Starting Separation", starting_separation,"xbar, rho = 1 , mu = ", mu, ", xbar = ", xbar)) + scale_y_log10()  + scale_x_log10() 


#print(Coalescence_Data_plot_ALL[,4])
}

print(p)

dev.off()
