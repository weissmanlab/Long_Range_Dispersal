library(ggplot2)
library(reshape2)


absorbing_boundary <- 999999

mutation_rate <- .01
ALPHA <- .5
SCALE_PARAMETER <- .005
#SCALE_PARAMETER <- 10
extra_scale <- 1/SCALE_PARAMETER #1/SCALE_PARAMETER
trial_cutoff <- 50
dist_1_num_of_trials_coalesced_or_absorbed <- 0
dist_100_num_of_trials_coalesced_or_absorbed <- 0
dist_10k_num_of_trials_coalesced_or_absorbed <- 0

dummy_df = read.table("sparse_pair_df.txt", sep = " ", header = TRUE)  # read text file 
dummy_df_coal_points = read.table("coalescence_points.txt", sep = " ", header = TRUE)  # read text file 
df <- data.frame(dummy_df)
df_coal_points <- data.frame(dummy_df_coal_points)

df$trial <- df$trial +1  # shift from c convention of starting from zero to r convention of starting from 1
df$time <- df$time +1 # shift starting generation from zero to one for log time axis
df$time <- mutation_rate*(df$time)
df$time <- df$time +1
df$time <- log10(df$time)


df$pos_1 <- asinh(df$pos_1*extra_scale*(mutation_rate)^(1/ALPHA))
df$pos_2 <- asinh(df$pos_2*extra_scale*(mutation_rate)^(1/ALPHA))

#print(df_coal_points$trial)
#if(length(df_coal_points$trial) > 0){
df_coal_points$trial <- df_coal_points$trial +1  # shift from c convention of starting from zero to r convention of starting from 1
df_coal_points$time <- df_coal_points$time +1 # shift starting generation from zero to one for log time axis
df_coal_points$time <- mutation_rate*(df_coal_points$time)
df_coal_points$time <- df_coal_points$time +1
df_coal_points$time <- log10(df_coal_points$time)

df_coal_points$pos_1 <- asinh(df_coal_points$pos_1*extra_scale*(mutation_rate)^(1/ALPHA))
df_coal_points$pos_2 <- asinh(df_coal_points$pos_2*extra_scale*(mutation_rate)^(1/ALPHA))




#}

#print(df_time_Test$time)

#y_max <- log10(max(df_time_Test$time))

df <- subset(df, alpha == ALPHA)
#print(nrow(df))
 df <- subset(df, scale_parameter == SCALE_PARAMETER)
 df <- subset(df, trial <= trial_cutoff)


df_coal_points <- subset(df_coal_points, alpha == ALPHA)
#print(nrow(df))
 df_coal_points <- subset(df_coal_points, scale_parameter == SCALE_PARAMETER)
 df_coal_points <- subset(df_coal_points, trial <= trial_cutoff)




df_time_Test <- df
df_time_Test <- subset(df_time_Test, pos_1 != "NA")
df_time_Test <- subset(df_time_Test, pos_2 != "NA")
#print(df_time_Test$pos_1)


y_max <- max(df_time_Test$time)

total_time <- max(df$time)
#print(total_time)

 df_distance_1 <- subset(df, distance == 1)

df_coal_points_distance_1 <- subset(df_coal_points, distance == 1)

 df_distance_10 <- subset(df, distance == 10)
  
 df_coal_points_distance_10 <- subset(df_coal_points, distance == 10)
 
 
df_distance_100 <- subset(df, distance == 100)
  
 df_coal_points_distance_100 <- subset(df_coal_points, distance == 100)

 df_distance_1000 <- subset(df, distance == 1000)
  
 df_coal_points_distance_1000 <- subset(df_coal_points, distance == 1000)

 
 df_distance_10k <- subset(df, distance == 10000)

df_coal_points_distance_10k <- subset(df_coal_points, distance == 10000)
#print(df_coal_points_distance_10k)



df_boundary = data.frame(right_boundary=rep(asinh(absorbing_boundary*extra_scale*(mutation_rate)^(1/ALPHA)), length(df$time)), left_boundary=rep(asinh(-1*absorbing_boundary*extra_scale*(mutation_rate)^(1/ALPHA)), length(df$time)), time = df$time )
#df_boundary = data.frame(right_boundary=rep(absorbing_boundary, length(df$time)), left_boundary=rep(-1*absorbing_boundary, length(df$time)), time = df$time )
#print(head(df_boundary ))

# if(ALPHA >= 2) {
#p <-  ggplot() + geom_path(data= df_distance_1, aes(pos_1, time, group = trial), color="RED", show.legend = FALSE) + geom_path(data= df_distance_1, aes(pos_2, time, group = trial), color="RED", show.legend = FALSE) + 
# geom_path(data= df_distance_100, aes(pos_1, time, group = trial), color="BLUE", show.legend = FALSE) + geom_path(data= df_distance_100, aes(pos_2, time, group = trial), color="BLUE", show.legend = FALSE) + 
# geom_path(data= df_distance_10k, aes(pos_1, time, group = trial), color="BLACK", show.legend = FALSE) + geom_path(data= df_distance_10k, aes(pos_2, time, group = trial), color="BLACK", show.legend = FALSE) + 
# theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank())   +  xlab("Arcsinh(x)") + ylab("Log(time)") 
#}

#if(ALPHA < 2 ) {
p <-  ggplot() + geom_path(data= df_distance_10, aes(pos_1, time, group = trial, alpha=I(.2)), color="RED", show.legend = FALSE) + geom_path(data= df_distance_10, aes(pos_2, time, group = trial, alpha=I(.2)), color="RED", show.legend = FALSE) +

 geom_path(data= df_distance_1000, aes(pos_1, time, group = trial), color="ORANGE", show.legend = FALSE) + geom_path(data= df_distance_1000, aes(pos_2, time, group = trial), color="ORANGE", show.legend = FALSE) + 
 geom_point(data= df_coal_points_distance_10, aes(pos_1, time), color="BLACK", show.legend = FALSE, shape = 4, size = 2) +
 geom_point(data= df_coal_points_distance_1000, aes(pos_1, time), color="blue1", show.legend = FALSE, shape = 4, size = 2) +
#geom_path(data= df_distance_10k, aes(pos_1, time, group = trial, alpha=I(0.2)), color="ORANGE", show.legend = FALSE) + geom_path(data= df_distance_10k, aes(pos_2, time, group = trial, alpha=I(0.2)), color="ORANGE", show.legend = FALSE) + 
#geom_point(data= df_coal_points_distance_10k, aes(pos_1, time), color="BLUE", show.legend = FALSE, shape = 4) +
theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + geom_path(data=df_boundary, aes(right_boundary, time), color='grey', show.legend = FALSE) + geom_path(data=df_boundary, aes(left_boundary, time), color = 'grey', show.legend = FALSE) + scale_y_continuous(limits = c(0, total_time))  + xlab("Arcsinh(x/xbar)") + ylab("Log(mutation*time)") + geom_hline(yintercept= 1, linetype="dashed", color = "black", size=.5)
 #+ ylab("Time") #+ ylab("Log(time)") 

#}

print(p)