#include <stdio.h>
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
{         if(argc != 6) {cout << "Wrong number of arguments.  Arguments are alpha, initial distance, number of trials, number of time steps, and timescale coarse graining." << endl; return 0;} 
  //coarse grain time... don't include every step - have another loop for coarse graining/smoothing time (10 steps for every recorded step)
  // longer timescales are necessary for free diffusion with no drift to create a proper distribution of coalescence times.  finit t_con seems to lead to a similar power law regardless of jump kernel- vary jump kernel and check this
  //int time = 0;
  // include periodic boundaries to ensure dist of coalescent times with jumps converges
//cout << fmod(-8, 3)<< endl;
  const int distance_off_set = 0;
  const int num_time_steps = atoi(argv[4]);
  const int time_scale_coarse_graining = atoi(argv[5]); //number of steps per (in between) recorded step -this is necessary so that we dont exceed memory requirements with arrays that are too large
  //const int num_mu_steps = 1000;  // number of mu increments in Laplace time 
  const int num_trials = atoi(argv[3]);
  const int num_distance_steps = 1;   //vary initial seperation exponentially for log plot of mean homozygosity as function of x for fixed mu
  const double periodic_boundary = 10000000; //position constrained between -pb and +pb
  const double D = 0;   // Diffusion constant
  const double cauchy_param = 0.00;//0.1;//0.01; //20; // controls cauchy power law jump kernel
  const double log_param = 0;   // controls lognormal jump kernel
  const double fisher_param = 1.00;
  const double cutoff = 0; // minimum jump size
  const double timestep = .1; // for deterministic drift term
  //const double mu_step = .0001; 
  const double t_con_inverse = .000;//.005; //.5 // (1/tcon) also for determinic drift term
  //const double rho_inverse = 10 ; // (1/rho) is for calculation of expectation over paths. Rho is population density.
  const double delta_function_width = 1;
  const double alpha = atof(argv[1]);;  // controls power law tail of jump kernel
  //long double position[num_time_steps];
  //long double Average_Position[num_time_steps][num_distance_steps] = {0};  // second index denotes initial position
  long double distance_list[num_distance_steps];
  //long double Contribution_from_each_trial[num_trials] = {0};
  //long double Contribution_from_each_trialEXPONENT[num_trials] = {0};
  //double dist_of_coalescent_times[num_time_steps] = {0};  //distribution of coalescent times for the given initial seperation 
  //double mean_homozygosity[num_mu_steps] = {0}; //probability of two individuals (lineages) being identical given initial seperation and mu
  
//long double dist_of_coalescent_times_ALL[num_time_steps][num_distance_steps] = {0};  //distribution of coalescent times for the given initial seperation 
  //long double mean_homozygosity_ALL[num_mu_steps][num_distance_steps] = {0}; //probability of two individuals (lineages) being identical given initial seperation and mu
  //long double normalization =0;
  //std::default_random_engine generator(time(0));
   std::mt19937 generator(time(0)); // mersenne twister psuedorandom number generator

  std::normal_distribution<double> norm_dist(0.0, 4*D);
  std::cauchy_distribution<double> cauchy_dist(0.0,1);
  std::lognormal_distribution<double> lognorm_dist(0.0,.5);
std::fisher_f_distribution<double> fisher_dist(2*alpha,2*alpha);
//use fisher distribution for tunable power law tail
  int dummy_counter = 0;

double initial_position = atof(argv[2]) ;  // initial signed distance between individuals
  double current_position;


          //char OUTPUTFILE2[50];
 // sprintf(OUTPUTFILE2, "alphavalue");
  std::stringstream file_name2;
         file_name2 <<  "alpha_value_"  <<  alpha  ;   // This is the directory name
         std::string stringfile2;
         file_name2 >> stringfile2; 
         
          std::stringstream file_name3;
         file_name3 <<  "mkdir alpha_value_"  << alpha ;   // This is mkdir directory name
         std::string stringfile3;
         getline(file_name3, stringfile3); 

char OUTPUTFILE99[50];
strcpy(OUTPUTFILE99, stringfile3.c_str());

cout << stringfile3 << endl;
system(OUTPUTFILE99);

char OUTPUTFILE100[50];
strcpy(OUTPUTFILE100, stringfile2.c_str());


chdir(OUTPUTFILE100);
for(int distance =0; distance < num_distance_steps; distance++)
{  //initial_position = 1; 
     initial_position = atof(argv[2]); 
   
      std::stringstream file_name4;
         file_name4 <<  "distance_value_"  <<  initial_position  ;   // This is the directory name
         std::string stringfile4;
         file_name4 >> stringfile4; 
         
          std::stringstream file_name5;
         file_name5 <<  "mkdir distance_value_"  << initial_position ;   // This is mkdir directory name
         std::string stringfile5;
         getline(file_name5, stringfile5); 

char OUTPUTFILE102[50];
strcpy(OUTPUTFILE102, stringfile5.c_str());

cout << stringfile5 << endl;
system(OUTPUTFILE102);

char OUTPUTFILE101[50];
strcpy(OUTPUTFILE101, stringfile4.c_str());


chdir(OUTPUTFILE101);

distance_list[distance] =  initial_position;
   current_position = fmod(initial_position, periodic_boundary);
for(int trial =0; trial < num_trials; trial++)
 {   //Contribution_from_each_trial[trial] = 0;
     //Contribution_from_each_trialEXPONENT[trial] = 0;
  

          char OUTPUTFILE[50];
  sprintf(OUTPUTFILE, "time_series");
  std::stringstream file_name;
         file_name <<  OUTPUTFILE  << "alpha" << alpha << "distance" << distance <<  "trial" << trial << ".txt" ;
         std::string stringfile;
         file_name >> stringfile; 
    ofstream fout7;
    fout7.open(stringfile);
  
   char OUTPUTFILE3[50];
  sprintf(OUTPUTFILE3, "entrance_and_exit_times");
  std::stringstream file_name3;
         file_name3 <<  OUTPUTFILE3  << "alpha" << alpha << "distance" << distance <<  "trial" << trial << ".txt" ;
         std::string stringfile3;
         file_name3 >> stringfile3; 
    ofstream fout8;
    fout8.open(stringfile3);
  bool inside_zone = false; 
   bool inside_zone_new; 
  //if(abs(current_position) <= delta_function_width)  // use step function with finite width as replacement for delta function
 //{inside_zone = true;}

  for (int time =0; time < (num_time_steps-1)*time_scale_coarse_graining; time++) {
    
     
 double signed_step_size =  sqrt(2*D)*norm_dist(generator) -timestep*t_con_inverse*current_position;
 
 double jump_size_cauchy = cauchy_dist(generator);
double jump_size_log = lognorm_dist(generator);
double jump_size_fisher = fisher_dist(generator);
//cauchy already takes both negative and psotive values //if(generator() > generator()) {jump_size_cauchy = - jump_size_cauchy;} // we want both positive and negative jumps
if(generator() > generator()) {jump_size_log = - jump_size_log;} // we want both positive and negative jumps
if(generator() > generator()) {jump_size_fisher = - jump_size_fisher;} // we want both positive and negative jumps
if(abs(jump_size_cauchy) > cutoff){signed_step_size = signed_step_size + cauchy_param*jump_size_cauchy;  } 
if(abs(jump_size_log) > cutoff){signed_step_size = signed_step_size + log_param*jump_size_log;  } 
if(abs(jump_size_fisher) > cutoff){signed_step_size = signed_step_size + fisher_param*jump_size_fisher;  }

//cout << jump_size_log << endl;
//cout << signed_step_size << endl;
//double dummy_time = atof(time);
//double time_index = int(floor(double(time)/time_scale_coarse_graining + .5));
// We coarse grain time by recording only every nth step.  This saves memory and allows us to exted to longer timescales.
//cout << int(floor(double(time)/time_scale_coarse_graining + .5))  << endl;
//dist_of_coalescent_times_ALL[int(floor(double(time)/time_scale_coarse_graining + .5))][distance] = dist_of_coalescent_times_ALL[int(floor(double(time)/time_scale_coarse_graining + .5))][distance] + Contribution_from_each_trial[trial]/num_trials;   // here we're adding up the contribution from each trial for a given time



 if(abs(current_position) <= delta_function_width)  // use step function with finite width as replacement for delta function
 { //Contribution_from_each_trial[trial] =  rho_inverse*exp(-Contribution_from_each_trialEXPONENT[trial]);
   //Contribution_from_each_trialEXPONENT[trial] += rho_inverse*timestep;  ;
    // Second term is the exponential discount factor which accounts for the probability that the two lineages have already coalesced.
    //fout8 << time << endl;
 

 inside_zone_new = true; 
 }

 else{ inside_zone_new = false; }

if(inside_zone_new == true && inside_zone == false)
{fout8 << time << " " ;}

if(inside_zone_new == false && inside_zone == true)
{fout8 << time << endl;}

inside_zone = inside_zone_new;


//cout << Contribution_from_each_trial[trial] << endl;
//cout << current_position << endl;
  fout7 << current_position << endl;
 //position[int(floor(double(time)/time_scale_coarse_graining + .5))] = current_position;
 //Average_Position[int(floor(double(time)/time_scale_coarse_graining + .5))][distance] += current_position/num_trials;
 current_position = fmod((current_position + signed_step_size),  periodic_boundary) ; 
    //cout << current_position << endl;
         //dummy_counter = dummy_counter +1;
        //cout << dummy_counter << endl;
       };
  fout7.close();
  fout8.close();
  current_position = fmod(initial_position, periodic_boundary);
  

  // next store position at every time step, plot as time series.
 

}

//normalize dist of coalescent times

/*
normalization = 0;

for (int time =0; time < num_time_steps; time++) {
  normalization = normalization + dist_of_coalescent_times_ALL[time][distance];
}

if(normalization != 0)
{for (int time =0; time < num_time_steps; time++) {
  dist_of_coalescent_times_ALL[time][distance] = dist_of_coalescent_times_ALL[time][distance]/normalization ;
   //dist_of_coalescent_times[time] =0;
}
normalization = 0;
}
else { cout << "normalization = 0 !" << endl;} 
*/
//Laplace transform dist of coalescent times to get mean homozygosity

/*
for (int mu =0; mu < num_mu_steps; mu++){
for (int time =0; time < num_time_steps; time++) {
    //mean_homozygosity_ALL[mu][distance] = mean_homozygosity_ALL[mu][distance] + dist_of_coalescent_times_ALL[time][distance]*exp(-mu*mu_step*time_scale_coarse_graining*time*timestep);
   mean_homozygosity_ALL[mu][distance] += dist_of_coalescent_times_ALL[time][distance]*exp(-mu*mu_step*time_scale_coarse_graining*time*timestep);
}}
*/


chdir("..");

}
chdir("..");
/*
ofstream fout;
fout.open("log_plot_mean_homozygosity_v_distance.txt");
for (int distance =0; distance < num_distance_steps; distance++) {
//mu = 10*mu_step
fout << log(mean_homozygosity_ALL[100][distance]) << endl;
  //fout << mean_homozygosity_ALL[100][distance] << endl;
  //fout <<  log(mean_homozygosity_ALL[999][distance]) << endl;

}
fout.close();
*/


/*
ofstream fout2;
fout2.open("log_plot_mean_homozygosity_v_distance.txt");
for (int time =0; time < num_time_steps; time++) {
fout2 << dist_of_coalescent_times_ALL[time][19] << endl;
 
}
fout2.close();
*/
/*for (int distance =0; distance < num_distance_steps; distance++) {

char OUTPUTFILE[50];
  sprintf(OUTPUTFILE, "time_series_averaged");
  std::stringstream file_name;
         file_name <<  OUTPUTFILE  << distance << ".txt" ;
         std::string stringfile;
         file_name >> stringfile; 
           

  ofstream fout4;
//fout4.open("time_series_averaged.txt");
fout4.open(stringfile);
for (int time =0; time < num_time_steps; time++) {
fout4 << Average_Position[time][distance] << endl;
 
}
fout4.close();
}

ofstream fout3;
fout3.open("distance_list.txt");
for (int distance =0; distance < num_distance_steps; distance++) {
fout3 <<  distance_list[distance] << endl;
 
}
fout3.close();
*/

/*
ofstream fout5;
fout5.open("Time_Series_Parameters.txt");

//fout5 <<  "num_distance_steps" << endl;
fout5 <<  num_distance_steps << endl;
//fout5 <<  "num_time_steps" << endl;
fout5 <<  num_time_steps << endl;
//fout5 <<  "time_scale_coarse_graining" << endl;
fout5 <<  time_scale_coarse_graining << endl;
//fout5 <<  "num_trials" << endl;
fout5 <<  num_trials << endl;
//fout5 <<  "timestep" << endl;
fout5 <<  timestep << endl;

//fout5 << "t_con_inverse"  << endl;
fout5 << t_con_inverse  << endl;

fout5.close();
*/
/*
const int distance_off_set = 0;
  const int num_time_steps = 1500;
  const int time_scale_coarse_graining = 10; //number of steps per (in between) recorded step -this is necessary so that we dont exceed memory requirements with arrays that are too large
  const int num_mu_steps = 1000;  // number of mu increments in Laplace time 
  const int num_trials = 50;
  const int num_distance_steps = 5;   //vary initial seperation exponentially for log plot of mean homozygosity as function of x for fixed mu
  const double periodic_boundary = 10000000; //position constrained between -pb and +pb
  const double D = 1;   // Diffusion constant
  const double cauchy_param = 0.00;//0.1;//0.01; //20; // controls cauchy power law jump kernel
  const double log_param = 0;   // controls lognormal jump kernel
  const double fisher_param = 0.00;
  const double cutoff = 0; // minimum jump size
  const double timestep = .1; // for deterministic drift term
  const double mu_step = .0001; 
  const double t_con_inverse = .000;//.005; //.5 // (1/tcon) also for determinic drift term
  const double rho_inverse = 10 ; // (1/rho) is for calculation of expectation over paths. Rho is population density.
  const double delta_function_width = 1;
  const double alpha = 0.0;  // controls power law tail of jump kernel
*/






/*
ofstream fout2;
fout2.open("dist_of_coalescent_times.txt");
for (int time =0; time < num_time_steps; time++) {
fout2 << dist_of_coalescent_times_ALL[time][num_distance_steps-1] << endl;
 
}
fout2.close();
*/
/*
ofstream fout3;
fout3.open("mean_homozygosity.txt");
for (int mu =0; mu < num_mu_steps; mu++) {
fout3 << mean_homozygosity_ALL[mu][14] << endl;
 
}
fout3.close();
*/

  return 0;

}
