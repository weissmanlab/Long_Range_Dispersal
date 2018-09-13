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
using namespace std;





int main(int argc, char* argv[])
{         if(argc != 6) {cout << "Wrong number of arguments.  Arguments are alpha, initial distance, number of trials, total number of time steps, rho_inverse " << endl; return 0;} 
  //coarse grain time... don't include every step - have another loop for coarse graining/smoothing time (10 steps for every recorded step)
  // longer timescales are necessary for free diffusion with no drift to create a proper distribution of coalescence times.  finit t_con seems to lead to a similar power law regardless of jump kernel- vary jump kernel and check this
  //int time = 0;
  // include periodic boundaries to ensure dist of coalescent times with jumps converges
//cout << fmod(-8, 3)<< endl;
  const int distance_off_set = 0;
  const int num_time_steps = atoi(argv[4]);
  const int time_scale_coarse_graining = 1;//atoi(argv[5]); //number of steps per (in between) recorded step -this is necessary so that we dont exceed memory requirements with arrays that are too large
  const int num_mu_steps = 10;  // number of mu increments in Laplace time 
  const int num_trials = atoi(argv[3]);
  const int num_distance_steps = 1;   //vary initial seperation exponentially for log plot of mean homozygosity as function of x for fixed mu
  const double periodic_boundary = 10000000; //position constrained between -pb and +pb
  const double D = 0;   // Diffusion constant
  //const double cauchy_param = 0.00;//0.1;//0.01; //20; // controls cauchy power law jump kernel
  //const double log_param = 0;   // controls lognormal jump kernel
  //const double fisher_param = 1.00;
  const double cutoff = 0; // minimum jump size
  const double timestep = 1; //const double timestep = .1; // for deterministic drift term and dist of coalescence.  for finite t_con this must be the same in part1 and part2
  const double mu_step = .01; 
  const double t_con_inverse = .000;//.005; //.5 // (1/tcon) also for determinic drift term
  const double rho_inverse = atof(argv[5]); // .1 ; // (1/rho) is for calculation of expectation over paths. Rho is population density.
  const double delta_function_width = 1;  //this must be same as in part 1
  const double alpha = atof(argv[1]);  // controls power law tail of jump kernel
  //long double position[num_time_steps];
  //long double Average_Position[num_time_steps][num_distance_steps] = {0};  // second index denotes initial position
  //long double distance_list[num_distance_steps];
  //long double Contribution_from_each_trial[num_trials] = {0};
  double *Contribution_from_each_trial = new double[num_trials];// {0};
  //long double Contribution_from_each_trialEXPONENT[num_trials] = {0};
  double *Contribution_from_each_trialEXPONENT= new double[num_trials];

for(int i =0; i < num_trials; i++)
{
   Contribution_from_each_trial[i] = 0;
   Contribution_from_each_trialEXPONENT[i] = 0;


}
  double *dist_of_coalescent_times = new double[num_time_steps];  //distribution of coalescent times for the given initial seperation 
 {for(int i =0; i < num_time_steps; i++)
{
   
   dist_of_coalescent_times[i] = 0;

}}


 double mean_homozygosity[num_mu_steps] = {0}; //probability of two individuals (lineages) being identical given initial seperation and mu
  
double mean_homozygosity_INDIVIDUAL_TRIAL[num_mu_steps] = {0};  

double mean_homozygosity_VARIANCE[num_mu_steps] = {0};  
//long double *dist_of_coalescent_times_ALL = new long double[num_time_steps][num_distance_steps];  //distribution of coalescent times for the given initial seperation 
  



   double mean_homozygosity_ALL[num_mu_steps][num_distance_steps] = {0}; //probability of two individuals (lineages) being identical given initial seperation and mu
  double normalization =0;
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

//cout << stringfile3 << endl;
//system(OUTPUTFILE99);

char OUTPUTFILE100[50];
strcpy(OUTPUTFILE100, stringfile2.c_str());


chdir(OUTPUTFILE100);
//for(int distance =0; distance < num_distance_steps; distance++)
//{  //initial_position = 1; 
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

//cout << stringfile5 << endl;
//system(OUTPUTFILE102);

char OUTPUTFILE101[50];
strcpy(OUTPUTFILE101, stringfile4.c_str());


chdir(OUTPUTFILE101);

//distance_list[distance] =  initial_position;
   current_position = fmod(initial_position, periodic_boundary);
for(int trial =0; trial < num_trials; trial++)
 {   //Contribution_from_each_trial[trial] = 0;
     //Contribution_from_each_trialEXPONENT[trial] = 0;
  for( int mu = 0; mu < num_mu_steps; mu++)
{
 mean_homozygosity_INDIVIDUAL_TRIAL[mu] =  0;


}
  

          char OUTPUTFILE88[50];
  //sprintf(OUTPUTFILE, "time_series");
  sprintf(OUTPUTFILE88, "entrance_and_exit_times");
  std::stringstream file_name88;
         file_name88 <<  OUTPUTFILE88  << "alpha" << alpha << "distance" << initial_position <<  "trial" << trial << ".txt" ;
         std::string stringfile88;
         file_name88 >> stringfile88; 
    ifstream fin88;
    //cout << stringfile88 << endl;
    fin88.open(stringfile88);
  //if(fin88.is_open() == false){cout << "NOT OPEN" << endl;}
  //if(fin88.is_open() == true){cout << "OPEN" << endl;}
int entrance_time = -1 ;  // If file is empty entrance and exit time will be the same and while loops will be ignored - dist of coalescent times will remain zero
  int exit_time= -1 ;
  //for (int time =0; time < (num_time_steps-1)*time_scale_coarse_graining; time++) {
    
    //for (int time =0; time < num_time_steps; time++) {
     for (int time =0; time < num_time_steps; time++) {



 

if (time < entrance_time)
{
   //dist_of_coalescent_times[time] += Contribution_from_each_trial[trial]/num_trials;
    Contribution_from_each_trialEXPONENT[trial] += 0;
    Contribution_from_each_trial[trial] = 0;
    dist_of_coalescent_times[time] += 0;  //no contribution from a trial when we are outside coalescence zone
}


if (time >= entrance_time && time < exit_time)
{
   Contribution_from_each_trialEXPONENT[trial] += rho_inverse*timestep;
//cout << Contribution_from_each_trial[trial] << endl;
Contribution_from_each_trial[trial] =  rho_inverse*exp(-Contribution_from_each_trialEXPONENT[trial]);
 dist_of_coalescent_times[time] += Contribution_from_each_trial[trial]/num_trials;

}

if (time >= exit_time)
{
   fin88 >> entrance_time >> exit_time;

}


       }
  



  fin88.close();
  
 

}











normalization = 0;

for (int time =0; time < num_time_steps; time++) {
  //normalization = normalization + dist_of_coalescent_times_ALL[time][distance];
   normalization +=  dist_of_coalescent_times[time];
}

if(normalization != 0)
{for (int time =0; time < num_time_steps; time++) {
  //dist_of_coalescent_times_ALL[time][distance] = dist_of_coalescent_times_ALL[time][distance]/normalization ;
   dist_of_coalescent_times[time] = dist_of_coalescent_times[time]/normalization ;

   //dist_of_coalescent_times[time] =0;
}
normalization = 0;
}
else { cout << "normalization = 0 !" << endl;} 

//Laplace transform dist of coalescent times to get mean homozygosity


for (int mu =0; mu < num_mu_steps; mu++){
for (int time =0; time < num_time_steps; time++) {
   
   mean_homozygosity[mu] += dist_of_coalescent_times[time]*exp(-mu*mu_step*time*timestep);

}}







// The loops above calculates the mean homozygosity

/*******************************************/
/*******************************************/

// The loops below calculates the variance and SD of the homozygosity














for(int trial =0; trial < num_trials; trial++)
 {   //Contribution_from_each_trial[trial] = 0;
     //Contribution_from_each_trialEXPONENT[trial] = 0;
  for( int mu = 0; mu < num_mu_steps; mu++)
{
 mean_homozygosity_INDIVIDUAL_TRIAL[mu] =  0;


}
  

          char OUTPUTFILE88[50];
  //sprintf(OUTPUTFILE, "time_series");
  sprintf(OUTPUTFILE88, "entrance_and_exit_times");
  std::stringstream file_name88;
         file_name88 <<  OUTPUTFILE88  << "alpha" << alpha << "distance" << initial_position <<  "trial" << trial << ".txt" ;
         std::string stringfile88;
         file_name88 >> stringfile88; 
    ifstream fin88;
    //cout << stringfile88 << endl;
    fin88.open(stringfile88);
  //if(fin88.is_open() == false){cout << "NOT OPEN" << endl;}
  //if(fin88.is_open() == true){cout << "OPEN" << endl;}
int entrance_time = -1 ;  // If file is empty entrance and exit time will be the same and while loops will be ignored - dist of coalescent times will remain zero
  int exit_time= -1 ;
  //for (int time =0; time < (num_time_steps-1)*time_scale_coarse_graining; time++) {
    
    //for (int time =0; time < num_time_steps; time++) {
     for (int time =0; time < num_time_steps; time++) {
 


 for( int mu = 0; mu < num_mu_steps; mu++)
{
 mean_homozygosity_INDIVIDUAL_TRIAL[mu] +=  Contribution_from_each_trial[trial]*exp(-mu*mu_step*time*timestep);


}




       }
  

for( int mu = 0; mu < num_mu_steps; mu++)
{
 
if(mean_homozygosity_INDIVIDUAL_TRIAL[0] != 0)
 {mean_homozygosity_INDIVIDUAL_TRIAL[mu] =  mean_homozygosity_INDIVIDUAL_TRIAL[mu]/mean_homozygosity_INDIVIDUAL_TRIAL[0];
  }

}


for( int mu = 0; mu < num_mu_steps; mu++)
{
 mean_homozygosity_VARIANCE[mu] +=  (mean_homozygosity_INDIVIDUAL_TRIAL[mu] - mean_homozygosity[mu])*(mean_homozygosity_INDIVIDUAL_TRIAL[mu] - mean_homozygosity[mu])/float(num_trials);
//cout << mean_homozygosity_INDIVIDUAL_TRIAL[mu] << endl;

}

  fin88.close();
  

}



/*******************************************/
/*******************************************/

//Now we output our results to files




char OUTPUTFILE[50];
  sprintf(OUTPUTFILE, "dist_of_coalescent_times_");
  std::stringstream file_name;
         file_name <<  OUTPUTFILE << "alpha_value_"<< alpha << "distance_value_" << initial_position << "rho_inverse_" << rho_inverse << ".txt" ;
         std::string stringfile;
         file_name >> stringfile; 
           

  ofstream fout4;
//fout4.open("time_series_averaged.txt");
fout4.open(stringfile);
for (int time =0; time < num_time_steps; time++) {
//fout4 << dist_of_coalescent_times_ALL[time][distance] << endl;
//if(time%500 ==0) {fout4 << time*timestep << " " << dist_of_coalescent_times[time] << endl;} // Coarse grain and only output every 500th value of the dist of coalescent times to save disk space (don't want a bunch of 5 million line files)

fout4 << time*timestep << " " << dist_of_coalescent_times[time] << endl; // Coarse grain and only output every 500th value of the dist of coalescent times to save disk space (don't want a bunch of 5 million line files)

//if(time !=0 && dist_of_coalescent_times[time] != dist_of_coalescent_times[time -1]) {fout4 << time*timestep << " " << dist_of_coalescent_times[time] << endl;}

}
fout4.close();


char OUTPUTFILE2[50];
  sprintf(OUTPUTFILE2, "mean_homozygosity_");
  std::stringstream file_name99;
         file_name99 <<  OUTPUTFILE2 << "alpha_value_"<< alpha << "distance_value_" << setw(7) << setfill('0') << initial_position  << "rho_inverse_" << rho_inverse << ".txt" ;
         std::string stringfile99;
         file_name99 >> stringfile99; 
           

  ofstream fout5;



fout5.open(stringfile99);
for (int mu =0; mu < num_mu_steps; mu++) {
fout5 << initial_position << " " << mu*mu_step << " " << mean_homozygosity[mu] << " " << (mean_homozygosity[mu] - sqrt(mean_homozygosity_VARIANCE[mu])/sqrt(float(num_trials))) <<  " " << (mean_homozygosity[mu] + sqrt(mean_homozygosity_VARIANCE[mu])/sqrt(float(num_trials))) << endl;
 
}
fout5.close();


chdir("..");

//}
chdir("..");


  return 0;

}
