#include <stdio.h> 
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>
#include <iostream>
#include <iomanip>
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
#include <type_traits>
#include <cstdlib>

using namespace std;

double periodic_boundary = 10000;  //range size is finite
double timestep = 1.0;
double delta_function_width = 1.0; // lattice spacing
double step_size_1; // used to increment random walk below
double step_size_2;
double step_x_1;
double step_y_1;
double step_x_2;
double step_y_2;
double normalized_step_x_1;
double normalized_step_y_1;
double normalized_step_x_2;
double normalized_step_y_2;
double mu_value;
double time_value;
double Contribution_from_each_trial;
double Contribution_from_each_trialEXPONENT;   
double prob_of_no_coalescence_yet; 



double return_rescaled_SD(double alpha, double scale_parameter)
{  double ORIG_SD  = alpha*sqrt( 2.0/( (alpha -2.0)*(alpha -1.0) )  ) ;
   return   sqrt(2)*sqrt(2)*(scale_parameter/ORIG_SD); // extra factor of sqrt(2) for multivariate normal
//ORIG_STD is the standard deviation of the original (two sided) F-distribution. 

// We rescale so that  sqrt(2)*sqrt(2)*scale_parameter is the standard deviation of the radial distribution we draw from. 
// This is needed for alpha > 2 only; in this regime dispersal is not a true Levy flight.  
// We probe these broader power laws with an F-distribution.
}
const double delta_function_area = 1.0; 
const double delta_function_height = 1.0; // this ensures coal kernel is normalized
const int num_mu_steps = 1;  // number of mu increments - should normally be 4
const double mu_step = 1.0; // should normally be .001 to recover all data 

 
double* return_1d_pointer(int index1_size )
{    double *POINTER_ARRAY_1D = new double[index1_size];
     for (int i = 0; i < index1_size; i++) {POINTER_ARRAY_1D[i] = 0;}    
     return POINTER_ARRAY_1D;
}



double** return_2d_pointer(int index1_size, int index2_size )
{    double **POINTER_ARRAY_2D = new double*[index1_size];
     for (int i = 0; i < index1_size; i++) {POINTER_ARRAY_2D[i] = new double[index2_size];}
     for(int i = 0; i < index1_size; i++){for(int j = 0; j < index2_size; j++){POINTER_ARRAY_2D[i][j] = 0; }}
     return POINTER_ARRAY_2D;
}
 
 
inline void print_Mean_Homozygosity(double alpha, double scale_parameter, double initial_position, int num_trials, int num_time_steps, double rho_inverse, double  mean_homozygosity[], double  lower_CI[], double upper_CI[])
{
char OUTPUTFILE2[50];
  sprintf(OUTPUTFILE2, "mean_homozygosity_2D_");
  std::stringstream file_name99;
         file_name99 <<  OUTPUTFILE2  << "scale_parameter_"<< scale_parameter << "_alpha_value_"<< alpha << "distance_value_" << setw(7) << setfill('0') << 
         initial_position  << "rho_inverse_" << rho_inverse << ".txt" ;
         std::string stringfile99;
         file_name99 >> stringfile99; 
  ofstream fout5;
fout5.open(stringfile99);
for (int mu =0; mu < num_mu_steps; mu++) {
fout5 << scale_parameter << " " << alpha <<  " " <<  rho_inverse  << " " << pow(10, mu)*mu_step << " " << initial_position << 
" " << mean_homozygosity[mu] << " " << lower_CI[mu] << " " <<  upper_CI[mu] << endl;
 // Here we output mean homozygosity as a function of mu and include error bars
}
fout5.close();


}
