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












