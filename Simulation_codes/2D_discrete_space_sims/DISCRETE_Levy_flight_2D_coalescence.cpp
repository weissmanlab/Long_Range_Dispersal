#include "DISCRETE_Levy_flight_2D_coalescence.h"
int main(int argc, char* argv[])
{         if(argc != 7) {cout << "Wrong number of arguments.  Arguments are alpha, initial distance, number of trials, total number of time steps, scale_parameter, rho_inverse " << endl; return 0;} 


const double alpha = atof(argv[1]);  // controls power law tail of jump kernel
const double initial_position = atof(argv[2]) ;  // initial signed distance between individuals
const    int num_trials = atoi(argv[3]); // number of trials
const    int num_time_steps = atoi(argv[4]); // number of time steps
const double  scale_parameter = atof(argv[5]); // Scale parameter, c,  of Levy stable jump kernel.  Note that c=(2* D_{\alpha}*timestep)^(1/alpha)
const double rescaled_SD =  return_rescaled_SD(alpha, scale_parameter); // This is for F -distribution when alpha > 2
const double rho_inverse = atof(argv[6]); 
const double prob_of_coalescence_this_timestep = 1- exp(-delta_function_height*rho_inverse*timestep);
double **List_of_single_trial_homozygosities = return_2d_pointer(num_mu_steps,  num_trials );
double mean_homozygosity[num_mu_steps] = {0}; //probability of two individuals being identical given initial seperation and mu


int num_bootstrap_samples = 1000; // DONT CHANGE WITHOUT CHANGING LINES 137-138 BELOW FOR PERCENTILES
double **List_of_bootstrap_mean_homozygosities = return_2d_pointer(num_mu_steps,  num_bootstrap_samples );


// set up random number generation for GSL
//////////////////////
//int seed = time(0); 
int seed = floor(pow(10,3)*initial_position*alpha); // deterministic seed that varies with input parameters
const gsl_rng_type * T;
T = gsl_rng_default;gsl_rng* r;
r = gsl_rng_alloc (T);
gsl_rng_set(r, seed); 
std::mt19937 generator(seed); 
/////////////////////

double current_position_x_1 = fmod(initial_position, periodic_boundary); ///pbc is at massive value (dbl_max).  It's just there to prevent numerical issues.
double current_position_y_1 = 0;

double current_position_x_2 = fmod(0.0, periodic_boundary); ///pbc is at massive value (dbl_max).  It's just there to prevent numerical issues.
double current_position_y_2 = 0;

current_position_x_1 = round(current_position_x_1);  // make sure positions are on the lattice
current_position_y_1 = round(current_position_y_1);
            
current_position_x_2 = round(current_position_x_2);
current_position_y_2 = round(current_position_y_2);


int separation_squared = 0;


for(int trial =0; trial < num_trials; trial++)
 { 
    Contribution_from_each_trial= 0;
    Contribution_from_each_trialEXPONENT= 0;   
    prob_of_no_coalescence_yet = 1; 

    for (int time =1; time < num_time_steps + 1; time++) 
    {  // Note that you must take at least one jump (and wait at least one generation) before you can coalesce
          
      time_value = double(time)*timestep;   
          
      if(alpha <= 2){
      
      step_size_1 =  gsl_ran_levy(r,  scale_parameter, alpha);
      gsl_ran_dir_2d_trig_method(r, &normalized_step_x_1, &normalized_step_y_1);
      //cout << normalized_step_x << " " <<  normalized_step_y << endl;
      step_x_1 = step_size_1*normalized_step_x_1;
      step_y_1 = step_size_1*normalized_step_y_1;
      
      step_size_2 =  gsl_ran_levy(r,  scale_parameter, alpha);
      gsl_ran_dir_2d_trig_method(r, &normalized_step_x_2, &normalized_step_y_2);
      //cout << normalized_step_x << " " <<  normalized_step_y << endl;
      step_x_2 = step_size_2*normalized_step_x_2;
      step_y_2 = step_size_2*normalized_step_y_2;
      
      }
      
       if(alpha == 2){
      step_x_1 = gsl_ran_levy(r,  scale_parameter, alpha);
      step_y_1 = gsl_ran_levy(r,  scale_parameter, alpha);
      
      
      step_x_2 = gsl_ran_levy(r,  scale_parameter, alpha);
      step_y_2 = gsl_ran_levy(r,  scale_parameter, alpha);
     
      }
      
      if(alpha > 2){ 
                     step_size_1 = rescaled_SD * gsl_ran_fdist(r, 2.0, 2*alpha );
                     gsl_ran_dir_2d_trig_method(r, &normalized_step_x_1, &normalized_step_y_1);
                     //cout << normalized_step_x << " " <<  normalized_step_y << endl;
                     step_x_1 = step_size_1*normalized_step_x_1;
                     step_y_1 = step_size_1*normalized_step_y_1;
                     
                     step_size_2 = rescaled_SD * gsl_ran_fdist(r, 2.0, 2*alpha );
                     gsl_ran_dir_2d_trig_method(r, &normalized_step_x_2, &normalized_step_y_2);
                     //cout << normalized_step_x << " " <<  normalized_step_y << endl;
                     step_x_2 = step_size_2*normalized_step_x_2;
                     step_y_2 = step_size_2*normalized_step_y_2;
        
        }
        //update position after latest step
            current_position_x_1 =  current_position_x_1 + step_x_1 ;   
            current_position_y_1 =  current_position_y_1 + step_y_1 ;   
            
            current_position_x_2 =  current_position_x_2 + step_x_2 ;   
            current_position_y_2 =  current_position_y_2 + step_y_2 ;   
            
            
            // apply periodic boundary conditions
            
            if(current_position_x_1 >= 0){ current_position_x_1 = fmod(current_position_x_1 + periodic_boundary/2.0, periodic_boundary) - periodic_boundary/2.0;  }
            if(current_position_x_1 < 0){ current_position_x_1 = fmod(current_position_x_1 - periodic_boundary/2.0, periodic_boundary) + periodic_boundary/2.0;  }
            
            if(current_position_y_1 >= 0){ current_position_y_1 = fmod(current_position_y_1 + periodic_boundary/2.0, periodic_boundary) - periodic_boundary/2.0;  }
            if(current_position_y_1 < 0){ current_position_y_1 = fmod(current_position_y_1 - periodic_boundary/2.0, periodic_boundary) + periodic_boundary/2.0;  }
            
             if(current_position_x_2 >= 0){ current_position_x_2 = fmod(current_position_x_2 + periodic_boundary/2.0, periodic_boundary) - periodic_boundary/2.0;  }
            if(current_position_x_2 < 0){ current_position_x_2 = fmod(current_position_x_2 - periodic_boundary/2.0, periodic_boundary) + periodic_boundary/2.0;  }
            
            if(current_position_y_2 >= 0){ current_position_y_2 = fmod(current_position_y_2 + periodic_boundary/2.0, periodic_boundary) - periodic_boundary/2.0;  }
            if(current_position_y_2 < 0){ current_position_y_2 = fmod(current_position_y_2 - periodic_boundary/2.0, periodic_boundary) + periodic_boundary/2.0;  }
            
        
            
            // round positions to nearest integers to keep trajectories on the lattice
            
            current_position_x_1 = round(current_position_x_1);
            current_position_y_1 = round(current_position_y_1);
            
            current_position_x_2 = round(current_position_x_2);
            current_position_y_2 = round(current_position_y_2);
            
            
            
            separation_squared = (current_position_x_2 - current_position_x_1)*(current_position_x_2 - current_position_x_1) + (current_position_y_2 - current_position_y_1)*(current_position_y_2 - current_position_y_1);
            
        //*******************************///////////////////// The chunk of code below checks if last step put lineages close enough to coalesce
        
        
        if( separation_squared == 0 )
        { 
        Contribution_from_each_trial = prob_of_no_coalescence_yet*prob_of_coalescence_this_timestep; //single trial dct
         prob_of_no_coalescence_yet *= exp(-delta_function_height*rho_inverse*timestep);  // update probability we havent coalesced
             for( int mu = 0; mu < num_mu_steps; mu++)
                { mu_value = pow(10, mu)*mu_step;
                  List_of_single_trial_homozygosities[mu][trial] +=         Contribution_from_each_trial*exp(-2.0*mu_value*time_value);

                } 
            
            
          }  

     };
  


    for( int mu = 0; mu < num_mu_steps; mu++)
                { 
                  mean_homozygosity[mu] += List_of_single_trial_homozygosities[mu][trial]/double(num_trials);
                } 
  
     
    current_position_x_1 = fmod(initial_position, periodic_boundary); //reset initial separation for each trial
    current_position_y_1 = 0;
    
    current_position_x_2 = fmod(0.0, periodic_boundary); //reset initial separation for each trial
    current_position_y_2 = 0;

 }
 

//next we use the poisson bootstrap for fast and efficient bootstrap sampling
for( int sample = 0; sample < num_bootstrap_samples; sample++){
    for(int trial =0; trial < num_trials; trial++){
        for( int mu = 0; mu < num_mu_steps; mu++){
                
                  List_of_bootstrap_mean_homozygosities[mu][sample] += gsl_ran_poisson( r, 1.0)*List_of_single_trial_homozygosities[mu][trial]/double(num_trials);
                  
                  
                }}}
                

// now sort the mean homozygosities from the bootstrap samples.  Bubble sort should be fine for 10^3 samples


double tmp1;
  double tmp2;
  for( int mu = 0; mu < num_mu_steps; mu++)
  { for(int sample1 = 0; sample1 < num_bootstrap_samples; sample1++)
    { for(int sample2 = sample1; sample2 < num_bootstrap_samples; sample2++)
      {    
        tmp1 = List_of_bootstrap_mean_homozygosities[mu][sample1]; 
        tmp2 = List_of_bootstrap_mean_homozygosities[mu][sample2];
        if( tmp2 < tmp1){  List_of_bootstrap_mean_homozygosities[mu][sample1] = tmp2; 
                           List_of_bootstrap_mean_homozygosities[mu][sample2] = tmp1;}
       }}
   }


                
//for( int sample = 1; sample < num_bootstrap_samples; sample++){

//cout << List_of_bootstrap_mean_homozygosities[0][sample] << endl;
//}
                
// percentile bootstrap
double lower_CI[num_mu_steps];
double upper_CI[num_mu_steps];

for( int mu = 0; mu < num_mu_steps; mu++){    
    lower_CI[mu] = List_of_bootstrap_mean_homozygosities[mu][159] ; //16th percentile
    upper_CI[mu] = List_of_bootstrap_mean_homozygosities[mu][839] ; //84th percentile
 
 } 


print_Mean_Homozygosity(alpha, scale_parameter, initial_position, num_trials, num_time_steps, rho_inverse, mean_homozygosity, lower_CI, upper_CI);

  return 0;
}
