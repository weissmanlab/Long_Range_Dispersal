#include "generate_sample_of_paths.h"
#include "assign_weights_to_paths.h"
#include "Calc_MH_simulations.h"
#include "Calc_MH_numerics.h"
#include "create_dataframe_for_plots.h"
#include "generate_plots.h"
#include <omp.h>
int main()
{ 
   



  for(int alpha_index = 0; alpha_index < 1; alpha_index++)
	{for(int mu_index = 1; mu_index < 2; mu_index++)
	 {

      //  #pragma omp parallel for schedule(auto) don't parallelize. it screws up the directory changes. Unavoidable race condition.   
	 	for(int distance_index = 0; distance_index < 14; distance_index++)
	      {  
             double ALPHA = 1.85;//1.25 + .2*double(alpha_index);
             double INIT_DISTANCE = exp(double(distance_index));
	      	 double MU = 1.0/pow(10.0, double(mu_index));
	      	 int NUM_TIME_STEPS = 250; //1e2; //int(10.0/MU);
             int NUM_TRIALS = 1e5 ;//int(double(1e4)/double(NUM_TIME_STEPS));
             double SCALE_PARAMETER = 250;
/*
             
             generate_sample_of_paths(ALPHA, INIT_DISTANCE, NUM_TRIALS, NUM_TIME_STEPS, SCALE_PARAMETER, MU);
             // generate_sample_of_paths(1.5, 1, 11, 11, 250, .1);
              assign_weights_to_paths(ALPHA, INIT_DISTANCE, NUM_TRIALS, NUM_TIME_STEPS, SCALE_PARAMETER, MU);
	      	
	      	for(int rho_inverse_index = 1; rho_inverse_index < 2; rho_inverse_index++)
	       	 { double RHO_INVERSE = 1.0/pow(10.0, double(rho_inverse_index));
              Calc_MH_simulations( ALPHA, INIT_DISTANCE, NUM_TRIALS, NUM_TIME_STEPS, MU, RHO_INVERSE);
               Calc_MH_numerics( ALPHA, INIT_DISTANCE,  NUM_TIME_STEPS,  SCALE_PARAMETER,  MU,  RHO_INVERSE);
         

	      	 
	       	 }

*/



	      }







	  }







	}

create_dataframe_for_plots();

for(int alpha_index = 0; alpha_index < 1; alpha_index++)
	{for(int mu_index = 1; mu_index < 2; mu_index++)
	      {  for(int rho_inverse_index = 1; rho_inverse_index < 2; rho_inverse_index++)
	       	 {  

	       	 	double ALPHA = 1.85; //1.25 + .2*double(alpha_index);
            double MU = 1.0/pow(10.0, double(mu_index));
	      	 double RHO_INVERSE = 1.0/pow(10.0, double(rho_inverse_index));  

                 generate_plots(ALPHA, MU, RHO_INVERSE);
	       	 
	       	 }}}




	return 0;
}