#include "part2.h"
int main(int argc, char* argv[])
{         if(argc != 7) {cout << "Wrong number of arguments.  Arguments are alpha, initial distance, number of trials, total number of time steps, scale_parameter, rho_inverse " << endl; return 0;} 
   
const double alpha = atof(argv[1]);  // controls power law tail of jump kernel
const double initial_position = atof(argv[2]) ;  // initial signed distance between individuals
const    int num_trials = atoi(argv[3]);
const    int num_time_steps = atoi(argv[4]);
const double  scale_parameter = atof(argv[5]); 
const double rho_inverse = atof(argv[6]); 
const double prob_of_coalescence_this_timestep = 1- exp(-delta_function_height*rho_inverse*timestep);
double *dist_of_coalescent_times = return_1d_pointer(num_trials);
double **List_of_single_trial_homozygosities = return_2d_pointer(num_mu_steps,  num_trials );
double **List_of_single_trial_CDF = return_2d_pointer(num_cdf_steps, num_trials );
  
chdir( alpha, initial_position, num_trials, num_time_steps, rho_inverse); // change directory
/*******************************/
//Now that we're in the proper directory with the entrance and exit times of interest, 
//we read them in and calculate the distribution of coalescence times and the mean homozygosity
/**************************************/
for(int trial =0; trial < num_trials; trial++)
 {  Contribution_from_each_trial= 0;
    Contribution_from_each_trialEXPONENT= 0;   
    prob_of_no_coalescence_yet = 1; 
          
    // Here we're opening the file containing entrance times associated with each individual trial         
    ifstream fin_entrance_and_exit_times;
    input_filename = return_input_filename(alpha, scale_parameter, initial_position, num_trials, num_time_steps, trial  );
    fin_entrance_and_exit_times.open(input_filename);

    /****************************************/ 
   
    double entrance_time = -1.0 ;  // If file is empty entrance and exit time will be the same and while loops will be ignored - dist of coalescent times will remain zero
    double exit_time= -1.0 ;

    for (int time =0; time < num_time_steps; time++) 
    {
        time_value = double(time)*timestep;  // physical time, rather than integer increments of the timestep

        if ( time_value >= entrance_time &&  time_value < exit_time) //update exponent and add to dist_of_coalescent_times[time] when in coalescence zone
        { 
         Contribution_from_each_trial = prob_of_no_coalescence_yet*prob_of_coalescence_this_timestep;
         dist_of_coalescent_times[time] += Contribution_from_each_trial/num_trials;  // average over trials as we go  
         prob_of_no_coalescence_yet *= exp(-delta_function_height*rho_inverse*timestep);  // update probability we havent coalesced
   
  
         for( int mu = 0; mu < num_mu_steps; mu++)
            { mu_value = pow(10, mu)*mu_step;
              List_of_single_trial_homozygosities[mu][trial] +=  Contribution_from_each_trial*exp(-2.0*mu_value*time_value);

            }
  

         for( int T = 0; T < num_cdf_steps; T++)
            { T_value = exp(double(T)/2.0)*timestep;      
              if( T_value > time_value ){ List_of_single_trial_CDF[T][trial] += Contribution_from_each_trial*timestep;}
             
            }

         }

           if (time_value >= exit_time){fin_entrance_and_exit_times >> entrance_time >> exit_time;} // read in new entrance and exit times

       }

     fin_entrance_and_exit_times.close();
     // bin histogram of single trial homozygosities here.
     Bin_single_trial_MH_and_CDF_into_histogram(num_trials, trial, List_of_single_trial_homozygosities, List_of_single_trial_CDF ); //, List_of_single_trial_homozygosities, List_of_single_trial_CDF);
}

LAPLACE_TRANSFORM_DCT_TO_MH(dist_of_coalescent_times, num_time_steps); // get mean homozygosity from dct
INTEGRATE_DCT_TO_CDF(dist_of_coalescent_times, num_time_steps); // get CDF from dct
/*******************************************/
// Calculate the confidence intervals using bootstrap
BOOTSTRAP_CI_MH_and_CDF(alpha, initial_position, num_trials, num_time_steps, rho_inverse, List_of_single_trial_homozygosities, List_of_single_trial_CDF);
//Now we output our results to files
print_Mean_Homozygosity(alpha, scale_parameter, initial_position, num_trials, num_time_steps, rho_inverse, mean_homozygosity, lower_CI, upper_CI);
print_CDF_of_coalescence_times(alpha, scale_parameter, initial_position, num_trials, num_time_steps, rho_inverse, CDF_of_DCT,  lower_CI_CDF, upper_CI_CDF);
print_histogram_of_single_trial_homozygosities(alpha, scale_parameter, initial_position, num_trials, num_time_steps, rho_inverse, hist_of_single_trial_homozygosities);
print_histogram_of_bootstrapped_mean_homozygosities(alpha, scale_parameter, initial_position, num_trials, num_time_steps, rho_inverse, hist_of_bootstrapped_mean_homozygosities);
print_Sorted_List_of_bootstrapped_mean_homozygosities(alpha, scale_parameter, initial_position, num_trials, num_time_steps, rho_inverse, Sorted_List_of_bootstrapped_mean_homozygosities);
print_list_of_single_trial_homozygosities(alpha, scale_parameter, initial_position, num_trials, num_time_steps, rho_inverse, List_of_single_trial_homozygosities);

chdir("..");
chdir("..");  // return to original directory
return 0;
}
