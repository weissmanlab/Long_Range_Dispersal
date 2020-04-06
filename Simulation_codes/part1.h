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

double periodic_boundary = DBL_MAX;  //massive value, range is effectively infinite.
double timestep = 1.0;
double delta_function_width = 1.0; //width of coalescence zone 
double signed_step_size; // used to increment random walk below

double return_rescaled_SD(double alpha, double scale_parameter)
{  double ORIG_SD  = 2*sqrt(alpha*(4*alpha -2)/((2*alpha -2)*(2*alpha -2)*(2*alpha -4)));
   return   (scale_parameter/ORIG_SD);
//ORIG_STD is the standard deviation of the original F-distribution. 

// We rescale so that scale_parameter is the standard deviation of the distribution we draw from. 
// This is needed for alpha > 2 only; in this regime dispersal is not a true Levy flight.  
// We probe these broader power laws with an F-distribution.
}


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



inline string return_output_filename(double alpha, double scale_parameter, double initial_position, int trial )
{ char INPUTFILE_entrance_and_exit_times[50];
  sprintf(INPUTFILE_entrance_and_exit_times, "entrance_and_exit_times");
  std::stringstream file_name_entrance_and_exit_times;
  file_name_entrance_and_exit_times <<  INPUTFILE_entrance_and_exit_times << "_scale_parameter_" << scale_parameter << "_alpha" << alpha << "distance" << initial_position <<  "trial" << trial << ".txt" ;
  std::string stringfile_entrance_and_exit_times;
  file_name_entrance_and_exit_times >> stringfile_entrance_and_exit_times; 
  return stringfile_entrance_and_exit_times;
}




