#include <stdio.h> 
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>
#include <iostream>
#include <unistd.h>
#include <cmath>
#include <random>
#include <iostream>
#include <string>
#include <string.h>
#include <fstream>
#include <sstream>
#include <stdlib.h>     /* atof */
#include <cfloat>
using namespace std;


 

 const double periodic_boundary = DBL_MAX;//100000000000; //position constrained between -pb and +pb 
  const double levy_param = 1.00;
  const double timestep = 1.0;//const double timestep = .1; // for deterministic drift term
  const double delta_function_width = 1.0; //atof(argv[5]);;
  // the scale parameter c is related to the generalized diffusion constant D as c =(4D*timestep)^(1/alpha)
  double jump_size_levy = 0;
  double jump_size_fisher = 0;

//ORIG_STD is the standard deviation of the fisher distribution.  We rescale so that scale_parameter is the standard deviation.




inline void chdir(double alpha, double initial_position, int num_trials, int num_time_steps, double scale_parameter)

{ std::stringstream file_name2;
         file_name2 <<  "alpha_value_"  <<  alpha  ;   // This is the directory name
         std::string stringfile2;
         file_name2 >> stringfile2;         
          std::stringstream file_name3;
         file_name3 <<  "mkdir alpha_value_"  << alpha ;   // This is mkdir directory name
         std::string stringfile3;
         getline(file_name3, stringfile3); 
char OUTPUTFILE_Create_Parent_Directory[50];
strcpy(OUTPUTFILE_Create_Parent_Directory, stringfile3.c_str());
system(OUTPUTFILE_Create_Parent_Directory);
char OUTPUTFILE_Parent_Directory[50];
strcpy(OUTPUTFILE_Parent_Directory, stringfile2.c_str());
chdir(OUTPUTFILE_Parent_Directory);
/********************************************/
   
      std::stringstream file_name4;
         file_name4 <<  "distance_value_"  <<  initial_position  ;   // This is the directory name
         std::string stringfile4;
         file_name4 >> stringfile4;   
          std::stringstream file_name5;
         file_name5 <<  "mkdir distance_value_"  << initial_position ;   // This is mkdir directory name
         std::string stringfile5;
         getline(file_name5, stringfile5); 
char OUTPUTFILE_Create_Child_Directory[50];
strcpy(OUTPUTFILE_Create_Child_Directory, stringfile5.c_str());
system(OUTPUTFILE_Create_Child_Directory);
char OUTPUTFILE_Child_Directory[50];
strcpy(OUTPUTFILE_Child_Directory, stringfile4.c_str());
chdir(OUTPUTFILE_Child_Directory);

}










