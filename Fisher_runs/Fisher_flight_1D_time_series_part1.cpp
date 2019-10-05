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
int main(int argc, char* argv[])
{         if(argc != 6) {cout << "Wrong number of arguments.  Arguments are alpha, initial distance, number of trials, total number of time steps, and scale parameter." << endl; return 0;}  //, and timescale coarse graining." << endl; return 0;} 
  // include periodic boundaries to ensure dist of coalescent times with jumps converges
 // const int distance_off_set = 0;
  const int num_time_steps = atoi(argv[4]);
 // const int time_scale_coarse_graining = 1; //atoi(argv[5]); //number of steps per (in between) recorded step -this is necessary so that we dont exceed memory requirements with arrays that are too large
  const int num_trials = atoi(argv[3]);
  //const int num_distance_steps = 1;   //vary initial seperation exponentially for log plot of mean homozygosity as function of x for fixed mu
  const double periodic_boundary = DBL_MAX;//100000000000; //position constrained between -pb and +pb
  //const double D = 0;   // Diffusion constant
 // const double cauchy_param = 0.00;//0.1;//0.01; //20; // controls cauchy power law jump kernel
 // const double log_param = 0;   // controls lognormal jump kernel
  const double levy_param = 1.00;
  //const double fisher_param = 0;//1.00;
  //const double cutoff = 0; // minimum jump size
  const double timestep = 1.0;//const double timestep = .1; // for deterministic drift term
  const double t_con_inverse = .000;//.005; //.5 // (1/tcon) also for determinic drift term
  const double delta_function_width = 1.0; //atof(argv[5]);;
  const double alpha = atof(argv[1]);  // controls power law tail of jump kernel
  const double scale_parameter = atof(argv[5]);//*pow(timestep, 1.0/alpha); // sets scale of levy alpha stable.  Coalescence zone "delta function" is of width one.  In order to test analytical predictions we want c >> 1.
  // the scale parameter c is related to the generalized diffusion constant D as c =(4D*timestep)^(1/alpha)
  
  double ORIG_STD = 2*sqrt(alpha*(4*alpha -2)/((2*alpha -2)*(2*alpha -2)*(2*alpha -4)));
  const gsl_rng_type * T;
  T = gsl_rng_default;
  gsl_rng* r;
   r = gsl_rng_alloc (T);
    gsl_rng_set(r, time(0));

  //std::default_random_engine generator(time(0));
   std::mt19937 generator(time(0)); // mersenne twister psuedorandom number generator
//std::normal_distribution<double> norm_dist(0.0, 4*D);
  //std::cauchy_distribution<double> cauchy_dist(0.0,1);
  //std::lognormal_distribution<double> lognorm_dist(0.0,.5);
//std::fisher_f_distribution<double> fisher_dist(2*alpha,2*alpha);
//std::Levy_alpha_stable_distribution<double> fisher_dist(2*alpha,2*alpha);
//double gsl_ran_levy(const gsl_rng *r, double c, double Alpha);

//double dummy;
//dummy = gsl_ran_levy(r,  1.5, 1.5);
//for (int i =0; i < 20; i++){cout << gsl_ran_levy(r,  1.5, 1.5)  << endl;}
//use fisher distribution for tunable power law tail
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
    fout_parameter_file << "number of trials " << num_trials << endl;
    fout_parameter_file << "number of time steps " << num_trials << endl;
    fout_parameter_file << "timestep size " << timestep << endl;
    fout_parameter_file << "t_con_inverse " << t_con_inverse << endl;
    fout_parameter_file << "scale_parameter " << scale_parameter << endl;
    fout_parameter_file.close();


for(int trial =0; trial < num_trials; trial++)
 {   //Contribution_from_each_trial[trial] = 0;
     //Contribution_from_each_trialEXPONENT[trial] = 0;
          char OUTPUTFILE[50];
  sprintf(OUTPUTFILE, "time_series");
  std::stringstream file_name;
         file_name <<  OUTPUTFILE  << "alpha" << alpha << "distance" << initial_position <<  "trial" << trial << ".txt" ;
         std::string stringfile;
         file_name >> stringfile; 
    ofstream fout7;
   // fout7.open(stringfile);
  
   char OUTPUTFILE3[50];
  sprintf(OUTPUTFILE3, "entrance_and_exit_times");
  std::stringstream file_name3;
         file_name3 <<  OUTPUTFILE3  << "alpha" << alpha << "distance" << initial_position <<  "trial" << trial << ".txt" ;
         std::string stringfile3;
         file_name3 >> stringfile3; 
    ofstream fout8;
    fout8.open(stringfile3);
  bool inside_zone = false; 
   bool inside_zone_new; 
  
  //for (int time =0; time < (num_time_steps-1)*time_scale_coarse_graining; time++) {
    
    for (int time =0; time < num_time_steps; time++) {
     
 double signed_step_size = 0;// sqrt(2*D)*norm_dist(generator) -timestep*t_con_inverse*current_position;
 
 //double jump_size_cauchy = cauchy_dist(generator);
//double jump_size_log = lognorm_dist(generator);
double jump_size_fisher = gsl_ran_fdist(r, 2*alpha, 2*alpha );
//double jump_size_levy = gsl_ran_levy(r,  scale_parameter, alpha);
//cauchy already takes both negative and psotive values //if(generator() > generator()) {jump_size_cauchy = - jump_size_cauchy;} // we want both positive and negative jumps
//if(generator() > generator()) {jump_size_log = - jump_size_log;} // we want both positive and negative jumps
if(generator() > generator()) {jump_size_fisher = - jump_size_fisher;} // we want both positive and negative jumps
//if(abs(jump_size_cauchy) > cutoff){signed_step_size = signed_step_size + cauchy_param*jump_size_cauchy;  } 
//if(abs(jump_size_log) > cutoff){signed_step_size = signed_step_size + log_param*jump_size_log;  } 
//if(abs(jump_size_fisher) > cutoff){signed_step_size +=  fisher_param*jump_size_fisher;  }



//cout << ORIG_STD << endl;
signed_step_size += (scale_parameter/ORIG_STD)*jump_size_fisher;

//ORIG_STD is the standard deviation of the fisher distribution.  We rescale so that scale_parameter is the standard deviation.
 
 if(abs(current_position) <= delta_function_width/2.0)  // use step function with finite width as replacement for delta function
 { inside_zone_new = true; }

 else{ inside_zone_new = false; }

if(inside_zone_new == true && inside_zone == false)
{fout8 << time << " " ;}

if(inside_zone_new == false && inside_zone == true)
{fout8 << time << endl;}

inside_zone = inside_zone_new;

 //if(time % 500 ==0){ 

 // fout7 << time*timestep << " " << current_position << endl;

//}
 
 current_position = fmod((current_position + signed_step_size),  periodic_boundary) ; 
   
       };
  //fout7.close();
  fout8.close();
  current_position = fmod(initial_position, periodic_boundary);
  

}

chdir("..");


chdir("..");

  return 0;

}
