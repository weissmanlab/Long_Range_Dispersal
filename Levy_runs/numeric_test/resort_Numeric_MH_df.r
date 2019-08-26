Numeric_Data_old <-   read.table("Numeric_MH_df.txt")
Numeric_Data_swapped <- Numeric_Data_old[, c(1,3,4,2,5,6,7)]
write.table(Numeric_Data_swapped,"Numeric_MH_df_swapped.txt",sep="\t",row.names=FALSE, col.names=FALSE)
#print(colnames(Numeric_Data_swapped) )
Numeric_Data_zeros <- Numeric_Data_swapped[Numeric_Data_swapped$V2 == 0,]
write.table(Numeric_Data_zeros,"numeric_MH_zeros.txt",sep="\t",row.names=FALSE,col.names=FALSE)
Numeric_Data_nonzero <- Numeric_Data_swapped[Numeric_Data_swapped$V2 > 0,]
Numeric_Data_nonzero[,4] <- log(Numeric_Data_nonzero[4])
write.table(Numeric_Data_nonzero,"numeric_MH_LOG_SCALE_ALL_but_zeros.txt",sep="\t",row.names=FALSE,col.names=FALSE)
