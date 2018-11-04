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
  
//We declare and initialize relevant variables in this section
/********************************/
 
  const int num_time_steps = atoi(argv[4]);
  const int num_mu_steps = 4;  // number of mu increments in Laplace time 
  const int num_trials = atoi(argv[3]);
  const int num_distance_steps = 1;   //vary initial seperation exponentially for log plot of mean homozygosity as function of x for fixed mu
  const double periodic_boundary = 10000000; //position constrained between -pb and +pb
  const double timestep = 1; //const double timestep = .1; // for deterministic drift term and dist of coalescence.  for finite t_con this must be the same in part1 and part2
  const double mu_step = .0001; 
  //const double t_con_inverse = .000;//.005; //.5 // (1/tcon) also for determinic drift term
  const double rho_inverse = atof(argv[5]); // .1 ; // (1/rho) is for calculation of expectation over paths. Rho is population density.
  const double delta_function_width = 1;  //this must be same as in part 1
  const double alpha = atof(argv[1]);  // controls power law tail of jump kernel
  double mean_homozygosity[num_mu_steps] = {0}; //probability of two individuals (lineages) being identical given initial seperation and mu
double mean_homozygosity_INDIVIDUAL_TRIAL[num_mu_steps] = {0}; // conditional expectation E[exp(-2mut)|path] ] = E[Hom|path]   
double second_moment_of_homozygosity_INDIVIDUAL_TRIAL[num_mu_steps] = {0}; // Second moment (conditional) E[exp(-4*mu*t)|path] ] = E[Hom^2|path]

double conditional_VARIANCE_of_homozygosity_INDIVIDUAL_TRIAL[num_mu_steps] = {0}; // conditional variance Var(exp(-2*mu*t)|path)  = E[Hom^2|path] - E[Hom|path]^2
double mean_homozygosity_VARIANCE_of_Conditional_Expectations[num_mu_steps] = {0};  // Var(E[Hom| path])
double mean_homozygosity_expectation_value_of_Conditional_VARIANCE[num_mu_steps] = {0};  // E[Var(Hom| path)]
double mean_homozygosity_TOTAL_VARIANCE[num_mu_steps] = {0} ;  // Var(exp(-2*mu*t)) = Var(Hom) = Var(E[Hom| path]) + E[Var(Hom| path)]
double initial_position = atof(argv[2]) ;  // initial signed distance between individuals
  double *Contribution_from_each_trial = new double[num_trials];// {0};
  double *Contribution_from_each_trialEXPONENT= new double[num_trials];
  double *dist_of_coalescent_times = new double[num_time_steps];
int Log_hist_of_single_trial_homozygosities[num_mu_steps][200] = {0};
for(int i =0; i < num_trials; i++)
{Contribution_from_each_trial[i] = 0; Contribution_from_each_trialEXPONENT[i] = 0;}
 for(int i =0; i < num_time_steps; i++)
{dist_of_coalescent_times[i] = 0;}
/********************************/
  //We declare and initialize relevant variables in the above section


//Next we change to the relevant directory
  
/**********************************/
  std::stringstream file_name2;
         file_name2 <<  "alpha_value_"  <<  alpha  ;   // This is the directory name
         std::string stringfile2;
         file_name2 >> stringfile2; 
         
       

char OUTPUTFILE_Parent_Directory[50];
strcpy(OUTPUTFILE_Parent_Directory, stringfile2.c_str());


chdir(OUTPUTFILE_Parent_Directory);

     initial_position = atof(argv[2]); 
   
      std::stringstream file_name4;
         file_name4 <<  "distance_value_"  <<  initial_position  ;   // This is the directory name
         std::string stringfile4;
         file_name4 >> stringfile4; 
         

char OUTPUTFILE_Child_Directory[50];
strcpy(OUTPUTFILE_Child_Directory, stringfile4.c_str());


chdir(OUTPUTFILE_Child_Directory);
/*******************************/
//Above we opened/moved into the relevant directory

//Now that we're in the proper directory with the entrance and exit times of interest, 
//we read them in and calculate the distribution of coalescence times and the mean homozygosity

/**************************************/
for(int trial =0; trial < num_trials; trial++)
 {  Contribution_from_each_trial[trial] = 0;
     Contribution_from_each_trialEXPONENT[trial] = 0;
          char INPUTFILE_entrance_and_exit_times[50];

  sprintf(INPUTFILE_entrance_and_exit_times, "entrance_and_exit_times");
  
  std::stringstream file_name_entrance_and_exit_times;
         file_name_entrance_and_exit_times <<  INPUTFILE_entrance_and_exit_times  << "alpha" << alpha << "distance" << initial_position <<  "trial" << trial << ".txt" ;
         std::string stringfile_entrance_and_exit_times;
         file_name_entrance_and_exit_times >> stringfile_entrance_and_exit_times; 
    ifstream fin_entrance_and_exit_times;
    /****************************************/
    // Here we're opening the file containing entrance times associated with each individual trial
    fin_entrance_and_exit_times.open(stringfile_entrance_and_exit_times);
 int entrance_time = -1 ;  // If file is empty entrance and exit time will be the same and while loops will be ignored - dist of coalescent times will remain zero
  int exit_time= -1 ;
  
for (int time =0; time < num_time_steps; time++) {
if (time >= entrance_time && time < exit_time) //update exponent and add to dist_of_coalescent_times[time] when in coalescence zone
{Contribution_from_each_trialEXPONENT[trial] += rho_inverse*timestep;
Contribution_from_each_trial[trial] =  rho_inverse*exp(-Contribution_from_each_trialEXPONENT[trial]);
 dist_of_coalescent_times[time] += Contribution_from_each_trial[trial]/num_trials;}

if (time >= exit_time)
{fin_entrance_and_exit_times >> entrance_time >> exit_time;}

}
  fin_entrance_and_exit_times.close();
}
//Laplace transform dist of coalescent times to get mean homozygosity
for (int mu =0; mu < num_mu_steps; mu++){
for (int time =0; time < num_time_steps; time++) {
   
   mean_homozygosity[mu] += dist_of_coalescent_times[time]*exp(-pow(10,mu)*2*mu_step*time*timestep); 
//factor of 2 in exponent is an important convention; we want the laplace transform L{P(t)}(2*mu) = E[exp(-2*mu*t)]
//For each trial we calculate the conditional expectation E[exp(-2*mu*t)|path] = E[Hom|path]
  //In order to correctly calculate the variance of the homozygosity we need  Var(exp(-2*mu*t)) = Var(Hom) = Var(E[Hom|path]) + E[Var(Hom|path)]  - Law of total variance

}}

// The loops above calculates the mean homozygosity

/*******************************************/
/*******************************************/

// The loops below calculates the variance and SD of the homozygosity

for(int trial =0; trial < num_trials; trial++)
 {   Contribution_from_each_trial[trial] = 0;
     Contribution_from_each_trialEXPONENT[trial] = 0;
       
       for( int mu = 0; mu < num_mu_steps; mu++)
     {
      mean_homozygosity_INDIVIDUAL_TRIAL[mu] =  0;
     second_moment_of_homozygosity_INDIVIDUAL_TRIAL[mu] = 0;
      conditional_VARIANCE_of_homozygosity_INDIVIDUAL_TRIAL[mu] = 0; 
     }
  

          char INPUTFILE_entrance_and_exit_times[50];
    // open the entrance and exit times for each trial again
  sprintf(INPUTFILE_entrance_and_exit_times, "entrance_and_exit_times");
  std::stringstream file_name_entrance_and_exit_times;
         file_name_entrance_and_exit_times <<  INPUTFILE_entrance_and_exit_times  << "alpha" << alpha << "distance" << initial_position <<  "trial" << trial << ".txt" ;
         std::string stringfile_entrance_and_exit_times;
         file_name_entrance_and_exit_times >> stringfile_entrance_and_exit_times; 
    ifstream fin_entrance_and_exit_times;
    fin_entrance_and_exit_times.open(stringfile_entrance_and_exit_times);
  int entrance_time = -1 ;  // If file is empty entrance and exit time will be the same and while loops will be ignored - dist of coalescent times will remain zero
  int exit_time= -1 ;
  
     for (int time =0; time < num_time_steps; time++) {
 
     if (time >= entrance_time && time < exit_time)
     {Contribution_from_each_trialEXPONENT[trial] += rho_inverse*timestep;
      Contribution_from_each_trial[trial] =  rho_inverse*exp(-Contribution_from_each_trialEXPONENT[trial]); 
//This array yields the conditional probability of coalescence for the current time, P(time|path).  By looping through num_timesteps we obtain and Laplace transform the distribution for all times
     // dist_of_coalescent_times[time] += Contribution_from_each_trial[trial]/num_trials;

      }

     if (time >= exit_time)
     {fin_entrance_and_exit_times >> entrance_time >> exit_time;
      }
          

// now we calculate the homozygosity for each trial
      for( int mu = 0; mu < num_mu_steps; mu++)
     {mean_homozygosity_INDIVIDUAL_TRIAL[mu] +=  Contribution_from_each_trial[trial]*exp(-pow(10, mu)*2*mu_step*time*timestep); // extra factor of 2 in exponent of Laplace transform is standard in pop gen
      
      second_moment_of_homozygosity_INDIVIDUAL_TRIAL[mu] += Contribution_from_each_trial[trial]*exp(-pow(10, mu)*4*mu_step*time*timestep); // E[Hom^2| path]
     }

      }
       

  // Then we determine the variance in the homozygosity
for( int mu = 0; mu < num_mu_steps; mu++)
{conditional_VARIANCE_of_homozygosity_INDIVIDUAL_TRIAL[mu] = second_moment_of_homozygosity_INDIVIDUAL_TRIAL[mu] - pow(mean_homozygosity_INDIVIDUAL_TRIAL[mu], 2);

mean_homozygosity_expectation_value_of_Conditional_VARIANCE[mu] += conditional_VARIANCE_of_homozygosity_INDIVIDUAL_TRIAL[mu]/num_trials;
mean_homozygosity_VARIANCE_of_Conditional_Expectations[mu] +=  pow((mean_homozygosity_INDIVIDUAL_TRIAL[mu] - mean_homozygosity[mu]),2)/float(num_trials); //
 




}
// Then we bin the homozygosities from each individual trial
// 
for( int mu = 0; mu < num_mu_steps; mu++)
{for(int QQ = 0; QQ < 200; QQ++)
   {
   if( mean_homozygosity_INDIVIDUAL_TRIAL[mu]  <= exp(-double(QQ)) && mean_homozygosity_INDIVIDUAL_TRIAL[mu] > exp(-double(QQ + 1))  )
       {  Log_hist_of_single_trial_homozygosities[mu][QQ] += 1 ;

       }
   

   }


}
 


 fin_entrance_and_exit_times.close();

}
/*******************************************/



//Now we output our results to files
/*******************************************/
char OUTPUTFILE[50];
  sprintf(OUTPUTFILE, "dist_of_coalescent_times_");
  std::stringstream file_name;
         file_name <<  OUTPUTFILE << "alpha_value_"<< alpha << "distance_value_" << initial_position << "rho_inverse_" << rho_inverse << ".txt" ;
         std::string stringfile;
         file_name >> stringfile; 
  ofstream fout4;

fout4.open(stringfile);
for (int time =0; time < num_time_steps; time++) {fout4 << time*timestep << " " << dist_of_coalescent_times[time] << endl;}
fout4.close();


char OUTPUTFILE2[50];
  sprintf(OUTPUTFILE2, "mean_homozygosity_");
  std::stringstream file_name99;
         file_name99 <<  OUTPUTFILE2 << "alpha_value_"<< alpha << "distance_value_" << setw(7) << setfill('0') << initial_position  << "rho_inverse_" << rho_inverse << ".txt" ;
         std::string stringfile99;
         file_name99 >> stringfile99; 

  ofstream fout5;

double SDOM = 0; // Standard deviation of the mean.  Defined in loop below.

fout5.open(stringfile99);
for (int mu =0; mu < num_mu_steps; mu++) {

mean_homozygosity_TOTAL_VARIANCE[mu] = mean_homozygosity_expectation_value_of_Conditional_VARIANCE[mu] + mean_homozygosity_VARIANCE_of_Conditional_Expectations[mu]; //law of total variance

SDOM = sqrt(mean_homozygosity_TOTAL_VARIANCE[mu])/sqrt(float(num_trials));

fout5 << initial_position << " " << pow(10, mu)*mu_step << " " << mean_homozygosity[mu] << " " << (mean_homozygosity[mu] - SDOM) <<  " " << (mean_homozygosity[mu] + SDOM) << endl;
 // Here we output mean homozygosity as a function of mu and include error bars - mean plus or minus standard deviation of the mean.
}
fout5.close();


char OUTPUTFILE3[50];
  sprintf(OUTPUTFILE3, "histogram_of_single_trial_homozygosities");
  std::stringstream file_name100;
         file_name100 <<  OUTPUTFILE3 << "alpha_value_"<< alpha << "distance_value_" << setw(7) << setfill('0') << initial_position  << "rho_inverse_" << rho_inverse << ".txt" ;
         std::string stringfile100;
         file_name100 >> stringfile100;

  ofstream fout6;

fout6.open(stringfile100);
for (int mu =0; mu < num_mu_steps; mu++) {
for (int QQ =0; QQ < 200; QQ++)
{fout6 << initial_position << " " << pow(10, mu)*mu_step << " " << -QQ  << " " << Log_hist_of_single_trial_homozygosities[mu][QQ] <<  endl;
 } // note that the log scale being used here is base e
// Here we output mean homozygosity as a function of mu and error bars - mean plus or minus standard deviation of the mean.
  }
  fout6.close();
 



chdir("..");
chdir("..");


return 0;
}
