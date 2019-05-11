#pragma once
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
#include<string.h>
#include <fstream>
#include <sstream>
#include <stdlib.h>     /* atof */
#include <type_traits>
#include <cstdlib>
#include <algorithm>
#include <vector>
//using namespace std;
void Calc_MH_simulations(const double ALPHA, const double INIT_DISTANCE, const int NUM_TRIALS, const int NUM_TIME_STEPS, const double MUTATION_RATE, const double RHO_INVERSE);
