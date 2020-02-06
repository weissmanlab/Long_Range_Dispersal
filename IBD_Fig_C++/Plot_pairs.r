library(ggplot2)
library(reshape2)


absorbing_boundary <- 999999

ALPHA <- 2
SCALE_PARAMETER <- 10

trial_cutoff <- 50
dist_1_num_of_trials_coalesced_or_absorbed <- 0
dist_100_num_of_trials_coalesced_or_absorbed <- 0
dist_10k_num_of_trials_coalesced_or_absorbed <- 0

dummy_df = read.table("sparse_pair_df.txt", sep = " ", header = TRUE)  # read text file 
df <- data.frame(dummy_df)


df$trial <- df$trial +1  # shift from c convention of starting from zero to r convention of starting from 1
df$time <- df$time +1 # shift starting generation from zero to one for log time axis
df$time <- log10(df$time)

df$pos_1 <- asinh(df$pos_1)
df$pos_2 <- asinh(df$pos_2)
#print(head(df))



#print(df_time_Test$time)

#y_max <- log10(max(df_time_Test$time))

df <- subset(df, alpha == ALPHA)
#print(nrow(df))
 df <- subset(df, scale_parameter == SCALE_PARAMETER)
 df <- subset(df, trial <= trial_cutoff)


df_time_Test <- df
df_time_Test <- subset(df_time_Test, pos_1 != "NA")
df_time_Test <- subset(df_time_Test, pos_2 != "NA")
#print(df_time_Test$pos_1)


y_max <- max(df_time_Test$time)

total_time <- max(df$time)
#print(total_time)

 df_distance_1 <- subset(df, distance == 1)



 df_distance_100 <- subset(df, distance == 100)
  
 
 
 
 
 df_distance_10k <- subset(df, distance == 10000)






df_boundary = data.frame(right_boundary=rep(asinh(absorbing_boundary), length(df$time)), left_boundary=rep(asinh(-1*absorbing_boundary), length(df$time)), time = df$time )

#print(head(df_boundary ))

# if(ALPHA >= 2) {
#p <-  ggplot() + geom_path(data= df_distance_1, aes(pos_1, time, group = trial), color="RED", show.legend = FALSE) + geom_path(data= df_distance_1, aes(pos_2, time, group = trial), color="RED", show.legend = FALSE) + 
# geom_path(data= df_distance_100, aes(pos_1, time, group = trial), color="BLUE", show.legend = FALSE) + geom_path(data= df_distance_100, aes(pos_2, time, group = trial), color="BLUE", show.legend = FALSE) + 
# geom_path(data= df_distance_10k, aes(pos_1, time, group = trial), color="BLACK", show.legend = FALSE) + geom_path(data= df_distance_10k, aes(pos_2, time, group = trial), color="BLACK", show.legend = FALSE) + 
# theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank())   +  xlab("Arcsinh(x)") + ylab("Log(time)") 
#}

#if(ALPHA < 2 ) {
p <-  ggplot() + geom_path(data= df_distance_1, aes(pos_1, time, group = trial, alpha=I(.2)), color="RED", show.legend = FALSE) + geom_path(data= df_distance_1, aes(pos_2, time, group = trial, alpha=I(.2)), color="RED", show.legend = FALSE) +
 #geom_path(data= df_distance_100, aes(pos_1, time, group = trial), color="BLUE", show.legend = FALSE) + geom_path(data= df_distance_100, aes(pos_2, time, group = trial), color="BLUE", show.legend = FALSE) + 
geom_path(data= df_distance_10k, aes(pos_1, time, group = trial, alpha=I(0.2)), color="BLACK", show.legend = FALSE) + geom_path(data= df_distance_10k, aes(pos_2, time, group = trial, alpha=I(0.2)), color="BLACK", show.legend = FALSE) + 
theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + geom_path(data=df_boundary, aes(right_boundary, time), color='grey', show.legend = FALSE) + geom_path(data=df_boundary, aes(left_boundary, time), color = 'grey', show.legend = FALSE) + scale_y_continuous(limits = c(0, y_max))  + xlab("Arcsinh(x)") + ylab("Log(time)")  #+ ylab("Time") #+ ylab("Log(time)") 

#}

print(p)