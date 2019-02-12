#DELETE OUTPUT FILE BEFORE REGENERATING.  
#IF THERE IS A TEXT FILE OTHER THAN THE INPUT FILES IN THE DIRECTORY YOU WILL GET AN ERROR

txt_files_ls = list.files( pattern="*.txt") 
# Read the files in, assuming comma separator
txt_files_df <- lapply(txt_files_ls, function(x) {read.table(file = x, header = F, sep =" ")})
# Combine them
combined_df <- do.call("rbind", lapply(txt_files_df, as.data.frame)) 
combined_df <- combined_df[order(combined_df$V1),]
rho_inverse <- .1
alpha <- 1.45
mu <- 0.001  # set whatever mu value you want
super_dummy = 0
sink(paste("outfile_ALL_alpha",alpha, "rho_inverse_", rho_inverse, sep=""))
for (i in 1:nrow(combined_df)){

if( combined_df[i, 2] == mu)
{   dummy <- paste0(alpha, " ", rho_inverse, " ", mu, " ", log(combined_df[i, 1]), " ", log(combined_df[i, 3]), " ", log(combined_df[i, 4]) , " ", log(combined_df[i, 5]) )
	#print(combined_df[i, 1], " ", combined_df[i, 3])
    #print(combined_df[i, 1])
     #print(combined_df[i, 3])
    cat(dummy)
     cat("\n")
     #plot( 1, 1)
     
     #print(super_dummy)
   #if(super_dummy != 0)
    # {points( combined_df[i, 1], combined_df[i, 3])
      #print("tester")
     #}


     #if(super_dummy == 0)
     #{plot( combined_df[i, 1], combined_df[i, 3])
      #super_dummy = 1
     #}

}



}


mu <- 0.01  # set whatever mu value you want

super_dummy = 0

for (i in 1:nrow(combined_df)){

if( combined_df[i, 2] == mu)
{   dummy <- paste0(alpha, " ", rho_inverse, " ", mu, " ", log(combined_df[i, 1]), " ", log(combined_df[i, 3]), " ", log(combined_df[i, 4]) , " ", log(combined_df[i, 5]) )
        #print(combined_df[i, 1], " ", combined_df[i, 3])
    #print(combined_df[i, 1])
     #print(combined_df[i, 3])
    cat(dummy)
     cat("\n")
     #plot( 1, 1)

     #print(super_dummy)
   #if(super_dummy != 0)
    # {points( combined_df[i, 1], combined_df[i, 3])
      #print("tester")
     #}


     #if(super_dummy == 0)
     #{plot( combined_df[i, 1], combined_df[i, 3])
      #super_dummy = 1
     #}

}



}




mu <- 0.1  # set whatever mu value you want

super_dummy = 0

for (i in 1:nrow(combined_df)){

if( combined_df[i, 2] == mu)
{   dummy <- paste0(alpha, " ", rho_inverse, " ", mu, " ", log(combined_df[i, 1]), " ", log(combined_df[i, 3]), " ", log(combined_df[i, 4]) , " ", log(combined_df[i, 5]) )
        #print(combined_df[i, 1], " ", combined_df[i, 3])
    #print(combined_df[i, 1])
     #print(combined_df[i, 3])
    cat(dummy)
     cat("\n")
     #plot( 1, 1)

     #print(super_dummy)
   #if(super_dummy != 0)
    # {points( combined_df[i, 1], combined_df[i, 3])
      #print("tester")
     #}


     #if(super_dummy == 0)
     #{plot( combined_df[i, 1], combined_df[i, 3])
      #super_dummy = 1
     #}

}



}





mu <- 1  # set whatever mu value you want

super_dummy = 0

for (i in 1:nrow(combined_df)){

if( combined_df[i, 2] == mu)
{   dummy <- paste0(alpha, " ", rho_inverse, " ", mu, " ", log(combined_df[i, 1]), " ", log(combined_df[i, 3]), " ", log(combined_df[i, 4]) , " ", log(combined_df[i, 5]) )
        #print(combined_df[i, 1], " ", combined_df[i, 3])
    #print(combined_df[i, 1])
     #print(combined_df[i, 3])
    cat(dummy)
     cat("\n")
     #plot( 1, 1)

     #print(super_dummy)
   #if(super_dummy != 0)
    # {points( combined_df[i, 1], combined_df[i, 3])
      #print("tester")
     #}


     #if(super_dummy == 0)
     #{plot( combined_df[i, 1], combined_df[i, 3])
      #super_dummy = 1
     #}

}



}






sink()
