#include "create_dataframe_for_plots.h"

void create_dataframe_for_plots()
{
        /*
    	//std::stringstream read_in_R_command;
         
        
         //read_in_R_command <<  "Rscript assign_weights_to_paths.r "  <<  ALPHA << " " << INIT_DISTANCE << " " <<  NUM_TRIALS << " " << NUM_TIME_STEPS << " " << SCALE_PARAMETER << " " <<  MUTATION_RATE;   
         // This produces the R command "Rscript assign_weights_to_paths.r ALPHA INIT_DISTANCE NUM_TRIALS NUM_TIME_STEPS SCALE_PARAMETER MUTATION_RATE"
         std::string dirname = "./Mh_plots";//read_in_R_command.str();
         
         char dirname_CHAR[200];   // conver to char array
		strcpy(dirname_CHAR, dirname.c_str());
         */


        chdir("./Mh_plots");
  		std::system("awk 1 mean_homozygosity_* > MH_data_frame.txt");

        chdir("..");

	
}