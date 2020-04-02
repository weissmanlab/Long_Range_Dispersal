library("minpack.lm")

library("abind")

library("ggplot2")

library("Rmisc")
library("grid")
library("dplyr")
library("libstableR")

set.seed(2)

#n <- 49000
n <- 46900
scale_parameter <- 1
time <- 1
#alpha = 1.45
#alpha = 2.0
alpha <- 1.6
pars <- c(alpha, 0, scale_parameter, 0)
rnd_x <- stable_rnd(n, pars)

x <- cumsum(sample(c(1, 1), n, TRUE))
#x <- cumsum(rnd_x)
rnd_y1 <- stable_rnd(n, pars)
y1 <- cumsum(rnd_y1) + 1000
set.seed(1)
rnd_y2 <- stable_rnd(n, pars)
y2 <- cumsum(rnd_y2) + 1500

set.seed(2)


#n <- 50500
alpha <- 2.0
pars <- c(alpha, 0, scale_parameter, 0)

rnd_y3 <- stable_rnd(n, pars)
y3 <- cumsum(rnd_y3)

set.seed(1)

rnd_y4 <- stable_rnd(n, pars)
y4 <- cumsum(rnd_y4) + 500

Y_sep = y2 - y1 
Time <- x
Space <- y1
df <- data.frame(Time,y1,y2, y3, y4)
ggplot(df, aes(Time)) +                    # basic graphical object
  geom_line(aes(y=Space), colour="red") +  # first layer
  geom_line(aes(y=y2), colour="red") +  # second layer
   geom_line(aes(y=y3), colour="orange") +
   geom_line(aes(y=y4), colour="orange") +
   #geom_segment(aes(x=0, xend=0, y=200, yend=0), arrow = arrow(length = unit(0.5, "cm"))) +
   #geom_segment(aes(x=0, xend=0, y=0, yend=200), arrow = arrow(length = unit(0.5, "cm"))) +
   geom_segment(aes(x=50000, xend=5000, y=-300, yend=-300), arrow = arrow(length = unit(0.5, "cm"))) +
   coord_flip() + 

   theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")
)



#plot(x, y1, type = "l")
#lines(x, y2)
