// In this code we sample trajectories that begin at a specified distance from the origin and end at the origin.  
//The number of timesteps for said trajectories is varied from one timestep to num_time_steps_max, 
//and for each number of timesteps a sample of num_trials trajectories is taken in order to calculate the distribution of coalescence times at that given time. 
// This procedure allows us to efficiently sample the paths that begin at an arbitrary point and end at the origin.
//#include <libstable/stable/src/stable.h> //  download this external library and add it to the your path (you should probably place it in the include directory.  Putting it in Long_Range_Dispersal is fine as well)
// MAKE SURE TO COMPILE libstable WITH make COMMAND - read the associated paper "libstable: Fast, Parallel, and High-Precision Computation of α-Stable Distributions in R, C/C++, and MATLAB"

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

using namespace std;
int main(int argc, char* argv[])
{      
    std::vector<double> v;
   if(argc != 6) {cout << "Wrong number of arguments.  Arguments are alpha, initial distance, number of trials for each terminal time  (we consider trajectories of varying time length) , maximum number of time steps, and scale parameter." << endl; return 0;}  //, and timescale coarse graining." << endl; return 0;} 
  // include periodic boundaries to ensure dist of coalescent times with jumps converges
 // const int distance_off_set = 0;
  const int num_time_steps_max = atoi(argv[4]);  //the total length of the trajectories considered
 // const int time_scale_coarse_graining = 1; //atoi(argv[5]); //number of steps per (in between) recorded step -this is necessary so that we dont exceed memory requirements with arrays that are too large
  const int num_trials = atoi(argv[3]); // numer of trials PER terminal timepoint. TOTAL
  const int Total_num_trials = num_time_steps_max*num_trials;
  //const int num_distance_steps = 1;   //vary initial seperation exponentially for log plot of mean homozygosity as function of x for fixed mu
  const double periodic_boundary = 10000000; //position constrained between -pb and +pb
  
  const double levy_param = 1.00;
  
  const double timestep = 1.0;//const double timestep = .1; // for deterministic drift term
  const double t_con_inverse = .000;//.005; //.5 // (1/tcon) also for determinic drift term
  const double delta_function_width = 1.0; //atof(argv[5]);;
  const double alpha = atof(argv[1]);  // controls power law tail of jump kernel
  const double scale_parameter = atof(argv[5]);//*pow(timestep, 1.0/alpha); // sets scale of levy alpha stable.  Coalescence zone "delta function" is of width one.  In order to test analytical predictions we want c >> 1.
  // the scale parameter c is related to the generalized diffusion constant D as c =(4D*timestep)^(1/alpha)
  const gsl_rng_type * T;
  T = gsl_rng_default;
  gsl_rng* r;
   r = gsl_rng_alloc (T);

 

   std::mt19937 generator(time(0)); // mersenne twister psuedorandom number generator

double gsl_ran_levy(const gsl_rng *r, double c, double Alpha);  // randomly generates numbers according to levy stable dist


  double initial_position = atof(argv[2]) ;  // initial signed distance between individuals
  double current_position;
//Relevant variables declared and initialized above
/********************************/
 
//Now we create and move into the directories in which we will output entrance and exit times
  std::stringstream file_name2;
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

     initial_position = atof(argv[2]); 
     //distance = initial_position;
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


   current_position = fmod(initial_position, periodic_boundary);
ofstream fout_parameter_file;
    fout_parameter_file.open("parameter_file.txt");
    fout_parameter_file << "alpha " << alpha << endl;
    fout_parameter_file << "distance " << initial_position << endl;
    fout_parameter_file << "number of trials per terminal timestep " << num_trials << endl;
    fout_parameter_file << "max number number of time steps considered " << num_time_steps_max << endl;
    fout_parameter_file << "max number number of time steps considered " << num_time_steps_max << endl;
    fout_parameter_file << "total number of trials " << Total_num_trials << endl;
    fout_parameter_file << "t_con_inverse " << t_con_inverse << endl;
    fout_parameter_file << "scale_parameter " << scale_parameter << endl;
    fout_parameter_file.close();

for(int num_time_steps =1; num_time_steps < num_time_steps_max; num_time_steps++)
{

for(int trial =0; trial < num_trials; trial++)
 {  
  
   char OUTPUTFILE3[50];
  sprintf(OUTPUTFILE3, "entrance_and_exit_times");
  std::stringstream file_name3;
         file_name3 <<  OUTPUTFILE3  << "alpha" << alpha << "distance" << initial_position << "num_time_steps" << num_time_steps <<  "trial" << trial << ".txt" ;
         std::string stringfile3;
         file_name3 >> stringfile3; 
    ofstream fout8;
    fout8.open(stringfile3);
 
  

 char OUTPUTFILE4[50];
  sprintf(OUTPUTFILE4, "free_and_contrained_CHOSEN_TIME_jump_sizes");
  std::stringstream file_name4;
         file_name4 <<  OUTPUTFILE4  << "alpha" << alpha << "distance" << initial_position << "num_time_steps" << num_time_steps <<  "trial" << trial << ".txt" ;
         std::string stringfile4;
         file_name4 >> stringfile4; 
    ofstream fout9;
    fout9.open(stringfile4);
  
  bool inside_zone = false; 
   bool inside_zone_new; 




    double free_trajectory[num_time_steps];  // dummy free trajectory generated will be used to construct constrained trajectory
    double constrained_trajectory[num_time_steps];  // the constrained trajectory we're interested in generating
    double free_step_list[num_time_steps];  // list of jumps taken at each timestep for free trajectory
    double constrained_step_list[num_time_steps]; // list of jumps taken at each timestep for constrained trajectory
 
    for (int time =0; time < num_time_steps; time++) {
     
 

double signed_step_size = gsl_ran_levy(r,  scale_parameter, alpha);
 
 free_step_list[time] = signed_step_size;  // list of jumps taken at each timestep for free trajectory
 free_trajectory[time] = current_position;
 constrained_step_list[time] = free_step_list[time]; // list of jumps taken at each timestep for constrained trajectory
 current_position = fmod((current_position + signed_step_size),  periodic_boundary) ; 
 
 //fout9 << free_step_list[time] << endl;
       }
  




  std::uniform_real_distribution<double> uniform_dist(0.0, num_time_steps);

  int chosen_time = int(floor(uniform_dist(generator)));  // Time at which we shift the trajectory by a large jump to ensure that the endpoint is the origin

  constrained_step_list[chosen_time] = free_step_list[chosen_time] - free_trajectory[num_time_steps -1]; // correct constrined list so that it will result in trajectory that ends at the origin

  fout9  << free_step_list[chosen_time] << " " << constrained_step_list[chosen_time] << endl;
  

 


  current_position = fmod(initial_position, periodic_boundary); //reset current position to inital position
     
    



     for (int time =0; time < num_time_steps; time++) {
     
 double signed_step_size = constrained_step_list[time];


 
 if(abs(current_position) <= delta_function_width/2.0)  // use step function with finite width as replacement for delta function
 { inside_zone_new = true; }

 else{ inside_zone_new = false; }

if(inside_zone_new == true && inside_zone == false)
{fout8 << time << " " ;}

if(inside_zone_new == false && inside_zone == true)
{fout8 << time << endl;}

inside_zone = inside_zone_new;

current_position = fmod((current_position + signed_step_size),  periodic_boundary) ;  
  

}

// The above process will generate trajectories that begin at the specified point and end at the origin.  However, the current sampling method is biased, 
//and we must weight/reweight the trajectories to correct for this bias.
// We want the probability of the constrained trajectory to be (proportional to) the probability of randomly drawing each of the constrained steps.
// Currently the probability of drawing a constrained trajectory is the probability of randomly drawing each of the free steps in the free step list.
// to correct for the bias we account for the one step that differs between the free and constrained step lists.
// We include with each sampled trajectory a weight equal to WEIGHT = Prob(constrained_step_list[chosen_time])/Prob(free_step_list[chosen_time])
// Including this weight in our calculations removes the bias and allows each path to contribute as if we had drawn directly from the conditional distribution.

// In part 1 we will output free_step_list[chosen_time] and constrained_step_list[chosen_time] along with the entrance and exit times
// In part 1.5 we will use and r code to evaluate the Levy stable dist at these two distance and use their ratrio to produce the weights


// In part two we will normalize the weights and perform a weighted average over paths to get the distribution of coalsecence times.  
// We then multiply the result by the probability of a trajectory ending at the origin to convert from the conditional expectation
// over constrained paths to the expectation over all paths.


fout8.close();
fout9.close();
current_position = fmod(initial_position, periodic_boundary); // reset to inital position

}}

chdir("..");


chdir("..");

  return 0;

}
