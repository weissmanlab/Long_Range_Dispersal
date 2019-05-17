This repository contains numeric and simulation data of coalescence for lineages undergoing long range dispersal.  

Data from the runs in the original format can be found in the directories new_Levy_runs and Fisher_runs (for Levy alpha Stable and Fisher distribution dispersal respectively).  These directories also contain c++ source codes and executables that produce the simulation data.  These codes work best with the g++ or clang compilers (other compilers have not been tested) and require the gsl library.  Examples of compilation and execution with the command line arguments used to produce our simulation results can be found in the example_Sim_Code_Runs.txt text file.

All of this data has been aggregated into dataframes.  One set of dataframes is for the mean homozygosity, the other set is for the distribution of coalescence times.

The first 3 columns of the mean homozygosity (MH) dataframes represent alpha, 1/rho, and mu respectively.  The fourth column is initial separation distance.  The fifth, sixth and seventh columns represent the MH and the lower and upper bounds of the associated bootstrap confidence intervals.  The MH values are stored on a linear scale in the simulation dataframes and on a natural log scale in the numeric data frames (except for zero initial separation - see details below).

The first two columns of the distribution of coalesene times (DCT) dataframes represent alpha and 1/rho.  The rest of the columns represent the DCT values at successive times, starting at t=0 and ending at t = 999.


Mean Homozygosity dataframes:
ALL_Fisher_runs_MH.txt -  contains all the MH data for the Fisher dist dispersal simulations only (stored on linear scale)
ALL_runs_MH.txt	-  ALL simulation MH data ( both distance and MH stored on linear scale)
numeric_MH_LOG_SCALE_ALL_but_zeros.txt - numeric values of MH for nonzero initial separation (both distance and MH stored on nautural log scale)
numeric_MH_zeros.txt - numeric values of MH for zero initial sepratation ( both distance and MH stored on linear scale)

DCT dataframes:
ALL_runs_DCT.txt - contains all the DCT simulation data (stored on linear scale for distance, time and DCT)
ALL_Fisher_runs_DCT.tx - contains all the DCT data for the Fisher dist dispersal simulations only (stored on linear scale for distance, time and DCT)



In the Plots directory you'll find subdirectories with R codes that produce specific plots found in the overleaf document. 
