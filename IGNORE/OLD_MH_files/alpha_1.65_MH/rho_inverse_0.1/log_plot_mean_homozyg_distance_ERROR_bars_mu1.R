library(Hmisc)

dummy <-   read.table("outfile_mu1")

#dist <- ts(dummy)
#plot.ts( dist)
#plot(dummy)
#abline(h = 4, v = 4, col = "gray60")
#abline(v = 4, col = "gray60")
#abline(a = 9.8, b = -2.95)
 distance <- as.vector(dummy[,1])
 class(distance)
 avg <- as.vector(dummy[,2])
 class(avg)
 upper <- dummy[,3]
 lower <- dummy[,4]
 plot(distance, avg, ylim=range(c(lower, upper)),
    pch=19)
arrows(distance, lower, distance, upper, length=0.05, angle=90, code=3)

 
 #plot(x, avg,
   # ylim=range(c(avg-sdev, avg+sdev)),
    #pch=19, xlab="Measurements", ylab="Mean +/- SD",
    #main="Scatter plot with std.dev error bars"
#)

 #with(dummy)
 #grid(nx=NA,ny=NULL)

