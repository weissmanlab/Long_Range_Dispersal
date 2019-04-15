// In this code we sample trajectories that begin at a specified distance from the origin and hit the origin at a random time.  

#pragma once

#include <stdio.h>
#include <vector>
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


void generate_sample_of_paths(const double ALPHA, const double INIT_DISTANCE, const int NUM_TRIALS, const int NUM_TIME_STEPS, const double SCALE_PARAMETER, const double MUTATION_RATE);
   
    
   