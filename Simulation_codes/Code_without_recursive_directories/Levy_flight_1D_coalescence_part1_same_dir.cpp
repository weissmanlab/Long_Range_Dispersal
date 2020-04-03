#include "part1_same_dir.h"
int main(int argc, char* argv[])
{         if(argc != 6) {cout << "Wrong number of arguments.  Arguments are alpha, initial distance, number of trials, total number of time steps, and scale parameter." << endl; return 0;} 

double alpha = atof(argv[1]);  // controls power law tail of jump kernel
double initial_position =atof(argv[2]);  // initial separation between lineages
int num_trials = atoi(argv[3]); 
int num_time_steps = atoi(argv[4]); 
double scale_parameter = atof(argv[5]); 
 // sets values of alpha, initial_position, num_trials, num_time_steps, scale_parameter
double ORIG_STD = 2*sqrt(alpha*(4*alpha -2)/((2*alpha -2)*(2*alpha -2)*(2*alpha -4))); // This is for F -distribution when alpha > 2

 int seed = time(0); const gsl_rng_type * T;T = gsl_rng_default;gsl_rng* r;r = gsl_rng_alloc (T);gsl_rng_set(r, time(0)); std::mt19937 generator(seed); 
//above line sets up random number generation

 // create and change into directory where we place the ouptut files for this simulation run
   double  current_position = fmod(initial_position, periodic_boundary);

for(int trial =0; trial < num_trials; trial++)
 {   

  //////////////////////////////////////////////  Create output file for this trial
  char OUTPUTFILE3[50];sprintf(OUTPUTFILE3, "entrance_and_exit_times");std::stringstream file_name3;
         file_name3 <<  OUTPUTFILE3  << "alpha" << alpha << "distance" << initial_position <<  "trial" << trial << ".txt" ;
         std::string stringfile3;file_name3 >> stringfile3; ofstream fout_trial_file;fout_trial_file.open(stringfile3);
////////////////////////////////

  bool inside_zone = false;  // this variable tells us whether or not we're close enough to coalesce
   bool inside_zone_new;  // this updates the above variable for the nex timestep

   for (int time =0; time < num_time_steps; time++) {  
         

          /////*******************************/////////////////////// The chunk of code below generates the next step in the random walk
    double signed_step_size = 0;

      if(alpha <= 2){jump_size_levy = gsl_ran_levy(r,  scale_parameter, alpha); signed_step_size +=  levy_param*jump_size_levy;}

      if(alpha > 2){ jump_size_fisher = gsl_ran_fdist(r, 2*alpha, 2*alpha );
        if(generator() > generator()) {jump_size_fisher = - jump_size_fisher;}  // F distribution is one sided.  We want both positive and negative jumps
      signed_step_size += (scale_parameter/ORIG_STD)*jump_size_fisher;
       }
       

        //*******************************///////////////////// The chunk of code below checks if last step put lineages close enough to coalesce

        if(abs(current_position) <= delta_function_width/2.0){ inside_zone_new = true; }  
        else{ inside_zone_new = false; }
        if(inside_zone_new == true && inside_zone == false){fout_trial_file << time << " " ;}
        if(inside_zone_new == false && inside_zone == true){fout_trial_file << time << endl;}
        inside_zone = inside_zone_new;
 
 /////////////////////////////////////////////////////////////////////////////////
          
          //update position after latest step
            current_position =  fmod((current_position + signed_step_size),  periodic_boundary) ;   // enforce pbc

         };
  
  fout_trial_file.close();  // close the file for this trial
  
current_position = fmod(initial_position, periodic_boundary);  //reset initial separation

}

  return 0;
}
