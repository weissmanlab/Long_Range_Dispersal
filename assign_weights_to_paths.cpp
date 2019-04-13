#include "assign_weights_to_paths.h"

void assign_weights_to_paths(const double ALPHA, const double INIT_DISTANCE, const int NUM_TRIALS, const int NUM_TIME_STEPS, const double SCALE_PARAMETER, const double MUTATION_RATE)
{ 

  

  std::stringstream read_in_R_command;
         
        
         read_in_R_command <<  "Rscript assign_weights_to_paths.r "  <<  ALPHA << " " << INIT_DISTANCE << " " <<  NUM_TRIALS << " " << NUM_TIME_STEPS << " " << SCALE_PARAMETER << " " <<  MUTATION_RATE;   
         // This produces the R command "Rscript assign_weights_to_paths.r ALPHA INIT_DISTANCE NUM_TRIALS NUM_TIME_STEPS SCALE_PARAMETER MUTATION_RATE"
         std::string R_command = read_in_R_command.str();
         //read_in_R_command >> R_command; 
         //std::cout << R_command << std::endl;
         char R_command_CHAR[200];   // conver to char array
strcpy(R_command_CHAR, R_command.c_str());
//std::cout << R_command_CHAR << std::endl;



  std::system(R_command_CHAR);

     





}