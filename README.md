This repository contains numeric and simulation data of coalescence for lineages undergoing Long range dispersal.  

Data from the runs in the original format can be found in the directories new_Levy_runs and Fisher runs (for Levy alpha Stable and Fisher distribution dispersal respectively).

This data has been aggregated into Dataframes.  One set of Dataframes is for the Mean homozygosity, the other is for the distribution of coalescence times.

The first 3 columns of the mean homozygosity dataframes represent alpha, 1/rho, and mu respectively.  The fourth column is initial separation distance.  The fifth sixth and seventh columns represent the Mean homozyogsity and the lower and upper bounds of the associated bootstrap confidence intervals.

The first two columns of the distribution of coalesene times (DCT) dataframes represent alpha and 1/rho.  The rest of the columns represent the DCT values at successive times, starting at t=0 and ending at t = 999.


In the Plots directory you'll find subdirectories with R codes that produce speficic plots found in the overleaf document. 
