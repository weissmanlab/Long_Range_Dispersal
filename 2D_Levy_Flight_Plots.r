library(ggplot2)
library(reshape2)
library(stabledist)
library(mvtnorm)
library(alphastable)
library(MASS)
library(wesanderson)




low_alpha_value <- .95
middle_alpha_value <- 1.5

mu <- .1
num_walks <- 1
num_gens <- 250
tail_cutoff <- 250

zero_matrix_all_alpha <- matrix(0, tail_cutoff, 3)

x_coord_array_all_alpha <- data.frame(zero_matrix_all_alpha)
y_coord_array_all_alpha <- data.frame(zero_matrix_all_alpha)


zero_matrix_single_timestep <- matrix(0, 1, 3)
x_coord_array_single_timestep <- data.frame(zero_matrix_single_timestep)
y_coord_array_single_timestep<- data.frame(zero_matrix_single_timestep)

x_coord_array_first_timestep <- data.frame(zero_matrix_single_timestep)
y_coord_array_first_timestep<- data.frame(zero_matrix_single_timestep)

x_coord_array_last_timestep <- data.frame(zero_matrix_single_timestep)
y_coord_array_last_timestep<- data.frame(zero_matrix_single_timestep)

for( QQ in 1:3)
{
	
if(QQ == 1){	
alpha <- low_alpha_value
scale_parameter <- 1*(mu)^(1/alpha)
}
if(QQ == 2){	
alpha <- middle_alpha_value 
scale_parameter <- 1*(mu)^(1/alpha)
}
if(QQ == 3){	
alpha <- 2
scale_parameter <- 1*(mu)^(1/alpha)
}
if(alpha > .5){
set.seed(seed = 135)
}

if(alpha == low_alpha_value){

#set.seed(seed = 138)
set.seed(seed = 141)
}
zero_matrix <- matrix(0, num_gens, 1)
x_coord_array <- data.frame(zero_matrix)
y_coord_array <- data.frame(zero_matrix)





rand_x_step_array <- data.frame(zero_matrix)
rand_y_step_array <- data.frame(zero_matrix)
if(alpha < 2){
dummy <- data.frame(mrstab.elliptical(num_gens,alpha,scale_parameter*matrix(c(1,.5,.5,1),2,2),c(0,0)))
}
if(alpha == 2){
dummy <- data.frame(rmvnorm(n=num_gens, mean=c(0,0), sigma =scale_parameter*matrix(c(1,.5,.5,1),2,2)))
}

rand_x_step_array <- dummy[,1]
rand_y_step_array <- dummy[,2]

x_coord_array <- cumsum(rand_x_step_array)
y_coord_array <- cumsum(rand_y_step_array)

x_coord_array_2 <- tail(x_coord_array, tail_cutoff)
y_coord_array_2 <- tail(y_coord_array, tail_cutoff)

constant_shift_x <- x_coord_array_2[1]
constant_shift_y <- y_coord_array_2[1]
if(QQ == 1)
{ x_coord_array_single_timestep[,1] <- x_coord_array_2[1/mu]  -  constant_shift_x 
  y_coord_array_single_timestep[,1] <- y_coord_array_2[1/mu]  -  constant_shift_y
  x_coord_array_first_timestep[,1] <- x_coord_array_2[1]  -  constant_shift_x 
  y_coord_array_first_timestep[,1] <- y_coord_array_2[1]  -  constant_shift_y
   x_coord_array_last_timestep[,1] <- x_coord_array_2[tail_cutoff]  -  constant_shift_x 
  y_coord_array_last_timestep[,1] <- y_coord_array_2[tail_cutoff]  -  constant_shift_y

  
  
for( TIME in 1:tail_cutoff)
{x_coord_array_2[TIME] <- x_coord_array_2[TIME]  -  constant_shift_x 
	y_coord_array_2[TIME] <- y_coord_array_2[TIME]  -  constant_shift_y
	 }
}

if(QQ == 2)
{x_coord_array_single_timestep[,2] <- x_coord_array_2[1/mu]  -  constant_shift_x -25 #-5 # - 7
  y_coord_array_single_timestep[,2] <- y_coord_array_2[1/mu]  -  constant_shift_y
   x_coord_array_first_timestep[,2] <- x_coord_array_2[1]  -  constant_shift_x -25
  y_coord_array_first_timestep[,2] <- y_coord_array_2[1]  -  constant_shift_y
   x_coord_array_last_timestep[,2] <- x_coord_array_2[tail_cutoff]  -  constant_shift_x  -25
  y_coord_array_last_timestep[,2] <- y_coord_array_2[tail_cutoff]  -  constant_shift_y
  
  
for( TIME in 1:tail_cutoff)
{x_coord_array_2[TIME] <- x_coord_array_2[TIME]  -  constant_shift_x  -25 #-5  #-7 #- 75
	y_coord_array_2[TIME] <- y_coord_array_2[TIME]  -  constant_shift_y
	 }
}


if(QQ == 3)
{x_coord_array_single_timestep[,3] <- x_coord_array_2[1/mu]  -  constant_shift_x -35#- 15
  y_coord_array_single_timestep[,3] <- y_coord_array_2[1/mu]  -  constant_shift_y
   x_coord_array_first_timestep[,3] <- x_coord_array_2[1]  -  constant_shift_x -35
  y_coord_array_first_timestep[,3] <- y_coord_array_2[1]  -  constant_shift_y
   x_coord_array_last_timestep[,3] <- x_coord_array_2[tail_cutoff]  -  constant_shift_x -35
  y_coord_array_last_timestep[,3] <- y_coord_array_2[tail_cutoff]  -  constant_shift_y
  
for( TIME in 1:tail_cutoff)
{x_coord_array_2[TIME] <- x_coord_array_2[TIME]  -  constant_shift_x -35 #- 15 #- 125
	y_coord_array_2[TIME] <- y_coord_array_2[TIME]  -  constant_shift_y 
	 }
}


if(alpha == low_alpha_value)
{ x_coord_array_all_alpha[,1] <- x_coord_array_2
  y_coord_array_all_alpha[,1] <- y_coord_array_2
	#print(x_coord_array_all_alpha[,1])
	
}

if(alpha == middle_alpha_value)
{  x_coord_array_all_alpha[,2] <- x_coord_array_2
   y_coord_array_all_alpha[,2] <- y_coord_array_2
	
	
}

if(alpha == 2)
{  x_coord_array_all_alpha[,3] <- x_coord_array_2
   y_coord_array_all_alpha[,3] <- y_coord_array_2
	
	
}

total_coord_array <- data.frame(matrix(0, num_gens, 2))
total_coord_array[,1] <- x_coord_array
total_coord_array[,2] <- y_coord_array


total_coord_array[,1] <- total_coord_array[,1] - total_coord_array[1,1]
total_coord_array[,2] <- total_coord_array[,2] - total_coord_array[1,2]


#print(total_coord_array)
#p <- p + geom_path(data=total_coord_array, aes(X1, X2), show.legend = FALSE) + geom_point(data=total_coord_array, aes(X1, X2), show.legend = FALSE) + theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
}
#print(p)

x_coord_array_all_alpha.long <- melt(x_coord_array_all_alpha)
y_coord_array_all_alpha.long <- melt(y_coord_array_all_alpha)

x_coord_array_single_timestep.long <- melt(x_coord_array_single_timestep)
y_coord_array_single_timestep.long <- melt(y_coord_array_single_timestep)
Total_coord_array_single_timestep <- data.frame(matrix(0, 3, 3))

Total_coord_array_single_timestep[,1] <- x_coord_array_single_timestep.long[,1]
Total_coord_array_single_timestep[,2] <- x_coord_array_single_timestep.long[,2]
Total_coord_array_single_timestep[,3] <- y_coord_array_single_timestep.long[,2]




x_coord_array_first_timestep.long <- melt(x_coord_array_first_timestep)
y_coord_array_first_timestep.long <- melt(y_coord_array_first_timestep)
Total_coord_array_first_timestep <- data.frame(matrix(0, 3, 3))

Total_coord_array_first_timestep[,1] <- x_coord_array_first_timestep.long[,1]
Total_coord_array_first_timestep[,2] <- x_coord_array_first_timestep.long[,2]
Total_coord_array_first_timestep[,3] <- y_coord_array_first_timestep.long[,2]


x_coord_array_last_timestep.long <- melt(x_coord_array_last_timestep)
y_coord_array_last_timestep.long <- melt(y_coord_array_last_timestep)
Total_coord_array_last_timestep <- data.frame(matrix(0, 3, 3))

Total_coord_array_last_timestep[,1] <- x_coord_array_last_timestep.long[,1]
Total_coord_array_last_timestep[,2] <- x_coord_array_last_timestep.long[,2]
Total_coord_array_last_timestep[,3] <- y_coord_array_last_timestep.long[,2]




#print(Total_coord_array_first_timestep)

Total_coord_array_all_alpha <- data.frame(matrix(0, 3*tail_cutoff, 3))
Total_coord_array_all_alpha[,1] <- x_coord_array_all_alpha.long[,1]
Total_coord_array_all_alpha[,2] <- x_coord_array_all_alpha.long[,2]
Total_coord_array_all_alpha[,3] <- y_coord_array_all_alpha.long[,2]
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#print(Total_coord_array_all_alpha)
p <- ggplot() +  geom_path(data=Total_coord_array_all_alpha, aes(X2, X3, color=X1), show.legend = FALSE) +  geom_point(data=Total_coord_array_single_timestep, aes(X2, X3), show.legend = FALSE, size = 2, shape = 17) +  geom_point(data=Total_coord_array_first_timestep, aes(X2, X3), show.legend = FALSE, size = 2) +  geom_point(data=Total_coord_array_last_timestep, aes(X2, X3), show.legend = FALSE, size = 2, shape = 15) + 
theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values=c("red", "blue", "orange"))
#scale_color_brewer(palette="Blues")




print(p)
