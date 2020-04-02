library(libstableR)
library(argparse)
library(mdatools)
library(ggplot2)
library(gridExtra)
library(Hmisc)
library(dplyr)

starting_separation <- 5
 distance <- 0
tbar <- 1



pdf(paste("Log_Plot_DCT_varying_alpha_short_distance_limit", "_tbar_val_", tbar,  ".pdf",sep=""))
p <- ggplot()
for(j in 1:3){
if(j==1)
{alpha <- .5}
if(j==2)
{alpha <- 1.5}
if(j==3)
{alpha <- 2}


#alpha <- 1.65




Final_Time <-1000


Numeric_approx_Data_plot_ALL <-  data.frame(matrix(0, Final_Time, 3))
for( i in 1:Final_Time )
{   Numeric_approx_Data_plot_ALL[i,1] <- i
	
	time_dummy <- i 
	 if(alpha == .5){
	 	dct <- gamma(1 + 1/alpha)*((time_dummy/tbar)^(-1/alpha))/(pi*tbar)
	 	
	 }
	 if(alpha == 1.5){
	 		dct <- (alpha -1)*((time_dummy/tbar)^(1/alpha -2))/(tbar/sin(pi/alpha))	 	
	 }
	 if(alpha == 2){
	 	dct <- (alpha -1)*((time_dummy/tbar)^(1/alpha -2))/(tbar/sin(pi/alpha))

	 	
	 }


	 
     
     
     
     	Numeric_approx_Data_plot_ALL[i,2] <-  dct
	Numeric_approx_Data_plot_ALL[i,3] <-  alpha
	
		}
Numeric_approx_Data_plot_ALL$X3 <- factor(Numeric_approx_Data_plot_ALL$X3)


colnames(Numeric_approx_Data_plot_ALL)[3] <- "Alpha" 



p <- p + geom_smooth(data=Numeric_approx_Data_plot_ALL, aes(x = X1, y =X2, color = Alpha, se = FALSE), linetype="dashed") + geom_point() + labs( x = "Time", y ="Dist of Coalescence Times p(t)") + ggtitle(paste("Starting Separation", 0,", rho = 1 , tbar = ", tbar)) + scale_y_log10()  + scale_x_log10() 


#print(Coalescence_Data_plot_ALL[,4])
}

print(p)

dev.off()
