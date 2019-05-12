#include "generate_plots.h"

void generate_plots(const double ALPHA, const double MU, const double RHO_INVERSE)
{




    std::stringstream read_in_R_command;
         
        
         read_in_R_command <<  "Rscript generate_plots.r "  <<  ALPHA << " " << MU << " " << RHO_INVERSE;   
         // This produces the R command "Rscript assign_weights_to_paths.r ALPHA INIT_DISTANCE NUM_TRIALS NUM_TIME_STEPS SCALE_PARAMETER MUTATION_RATE"
         std::string R_command = read_in_R_command.str();
         
         char R_command_CHAR[200];   // conver to char array
    strcpy(R_command_CHAR, R_command.c_str());




     std::system(R_command_CHAR);

     








}