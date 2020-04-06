#include "part1_same_dir.h"
int main(int argc, char* argv[])
{ if(argc != 6) {cout << "Wrong number of arguments.  Arguments are alpha, initial distance, number of trials, total number of time steps, and scale parameter." << endl; return 0;} 
const double alpha = atof(argv[1]);  // controls power law tail of jump kernel
const double initial_position =atof(argv[2]);  // initial separation between lineages
const    int num_trials = atoi(argv[3]);  // number of trials
const    int num_time_steps = atoi(argv[4]);  // number of time steps
const double scale_parameter = atof(argv[5]); // Scale parameter, c,  of Levy stable jump kernel.  Note that =(2* D_{\alpha}*timestep)^(1/alpha)
const double rescaled_SD =  return_rescaled_SD(alpha, scale_parameter); // This is for F -distribution when alpha > 2

// set up random number generation for GSL
//////////////////////
int seed = time(0); 
const gsl_rng_type * T;
T = gsl_rng_default;gsl_rng* r;
r = gsl_rng_alloc (T);
gsl_rng_set(r, seed); 
std::mt19937 generator(seed); 
/////////////////////

//chdir(alpha, initial_position, num_trials, num_time_steps, scale_parameter); // create and change into directory where we want output files
double  current_position = fmod(initial_position, periodic_boundary); //pbc is at massive value (dbl_max).  It's just there to prevent numerical issues.

for(int trial =0; trial < num_trials; trial++)
 {  ofstream fout_trial_file;
    string stringfile = return_output_filename( alpha, scale_parameter, initial_position,  trial);
    fout_trial_file.open(stringfile);
    bool inside_zone = false;  // this variable tells us whether or not we're close enough to coalesce during a given timestep
    bool inside_zone_new;  // this updates the above variable for the nex timestep

    for (int time =0; time < num_time_steps; time++) {  
          
      if(alpha <= 2){ signed_step_size =  gsl_ran_levy(r,  scale_parameter, alpha);}
      
      if(alpha > 2){ signed_step_size = rescaled_SD * gsl_ran_fdist(r, 2*alpha, 2*alpha );
        if(gsl_ran_bernoulli( r, .5) > 0) {signed_step_size = - signed_step_size;}  // F distribution is one sided.  We want two sided jumps.
        }
       
        //*******************************///////////////////// The chunk of code below checks if last step put lineages close enough to coalesce
        if(abs(current_position) <= delta_function_width/2.0){ inside_zone_new = true; }  
        else{ inside_zone_new = false; }
        
        if(inside_zone_new == true && inside_zone == false){fout_trial_file << time << " " ;}  // this is an entrance time
        if(inside_zone_new == false && inside_zone == true){fout_trial_file << time << endl;}  // this is an exit time
        
        inside_zone = inside_zone_new;  // update for next time step
 
 
          //update position after latest step
            current_position =  fmod((current_position + signed_step_size),  periodic_boundary) ;   // enforce pbc

         };
  
    fout_trial_file.close();  // close the file for this trial
  
    current_position = fmod(initial_position, periodic_boundary);  //reset initial separation for each trial

   }
chdir("..");
chdir("..");  // return to original directory
  return 0;
}
