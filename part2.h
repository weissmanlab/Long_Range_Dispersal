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


 const double delta_function_width = 1.0; //atof(argv[5]);  //this must be same as in part 1
 const double delta_function_height = 1.0/delta_function_width; 
const int num_mu_steps = 7;  // number of mu increments in Laplace time 
const int num_cdf_steps = 29;
const int num_distance_steps = 1;   //vary initial seperation exponentially for log plot of mean homozygosity as function of x for fixed mu
const double timestep = 1.0; //const double timestep = .1; // for deterministic drift term and dist of coalescence.  for finite t_con this must be the same in part1 and part2
const double mu_step = .001;  // .0001; //change back later 
const int num_histogram_bins = 2000;
double hist_of_single_trial_homozygosities[num_mu_steps][num_histogram_bins] = {0};
double hist_of_bootstrapped_mean_homozygosities[num_mu_steps][num_histogram_bins] = {0};
const int num_samples_bootstrapped_means = 10000;
double List_of_bootstrapped_mean_homozygosities[num_mu_steps][num_samples_bootstrapped_means] = {0};
double Sorted_List_of_bootstrapped_mean_homozygosities[num_mu_steps][num_samples_bootstrapped_means] = {0};


double hist_of_single_trial_cdfs[num_cdf_steps][num_histogram_bins] = {0};
double hist_of_bootstrapped_cdfs[num_cdf_steps][num_histogram_bins] = {0};
const int num_samples_bootstrapped_cdfs = 10000;
double List_of_bootstrapped_cdfs[num_cdf_steps][num_samples_bootstrapped_cdfs] = {0};
double Sorted_List_of_bootstrapped_cdfs[num_cdf_steps][num_samples_bootstrapped_cdfs] = {0};

inline void chdir(double alpha, double initial_position, int num_trials, int num_time_steps, double rho_inverse)

{std::stringstream file_name2;
         file_name2 <<  "alpha_value_"  <<  alpha  ;   // This is the directory name
         std::string stringfile2;
         file_name2 >> stringfile2; 
char OUTPUTFILE_Parent_Directory[50];
strcpy(OUTPUTFILE_Parent_Directory, stringfile2.c_str());
chdir(OUTPUTFILE_Parent_Directory);
      std::stringstream file_name4;
         file_name4 <<  "distance_value_"  <<  initial_position  ;   // This is the directory name
         std::string stringfile4;
         file_name4 >> stringfile4; 
char OUTPUTFILE_Child_Directory[50];
strcpy(OUTPUTFILE_Child_Directory, stringfile4.c_str());
chdir(OUTPUTFILE_Child_Directory);

}






inline void print_list_of_single_trial_homozygosities(double alpha, double initial_position, int num_trials, int num_time_steps, double rho_inverse, double** List_of_single_trial_homozygosities  )
{
char OUTPUTFILE6[50];
  sprintf(OUTPUTFILE6, "list_of_single_trial_homozygosities");
  std::stringstream file_name103;
         file_name103 <<  OUTPUTFILE6 << "alpha_value_"<< alpha << "distance_value_" << setw(7) << setfill('0') << initial_position  << "rho_inverse_" << rho_inverse << ".txt" ;
         std::string stringfile103;
         file_name103 >> stringfile103;

  ofstream fout9;

fout9.open(stringfile103);


for (int mu =0; mu < num_mu_steps; mu++) {
for (int QQ =0; QQ < num_trials; QQ++)
{fout9 << initial_position << " " << pow(10, mu)*mu_step << " " << List_of_single_trial_homozygosities[mu][QQ]  << " " << 1.0/double(num_trials)<<  endl;
 } 
// Here we output mean homozygosity as a function of mu and error bars - mean plus or minus standard deviation of the mean.
  }
  fout9.close();
}








inline void print_Sorted_List_of_bootstrapped_mean_homozygosities(double alpha, double initial_position, int num_trials, int num_time_steps, double rho_inverse, double Sorted_List_of_bootstrapped_mean_homozygosities[][num_samples_bootstrapped_means]  )
{
char OUTPUTFILE5[50];
  sprintf(OUTPUTFILE5, "sorted_list_of_bootstrapped_mean_homozygosities");
  std::stringstream file_name102;
         file_name102 <<  OUTPUTFILE5 << "alpha_value_"<< alpha << "distance_value_" << setw(7) << setfill('0') << initial_position  << "rho_inverse_" << rho_inverse << ".txt" ;
         std::string stringfile102;
         file_name102 >> stringfile102;

  ofstream fout8;

fout8.open(stringfile102);


for (int mu =0; mu < num_mu_steps; mu++) {
for (int QQ =0; QQ < num_samples_bootstrapped_means; QQ++)
{fout8 << initial_position << " " << pow(10, mu)*mu_step << " " << Sorted_List_of_bootstrapped_mean_homozygosities[mu][QQ]  << " " << 1.0/double(num_samples_bootstrapped_means)<<  endl;
 } 
// Here we output mean homozygosity as a function of mu and error bars - mean plus or minus standard deviation of the mean.
  }
  fout8.close();




}






inline void print_histogram_of_bootstrapped_mean_homozygosities(double alpha, double initial_position, int num_trials, int num_time_steps, double rho_inverse, double hist_of_bootstrapped_mean_homozygosities[][num_histogram_bins]  )
{
char OUTPUTFILE4[50];
  sprintf(OUTPUTFILE4, "histogram_of_bootstrapped_mean_homozygosities");
  std::stringstream file_name101;
         file_name101 <<  OUTPUTFILE4 << "alpha_value_"<< alpha << "distance_value_" << setw(7) << setfill('0') << initial_position  << "rho_inverse_" << rho_inverse << ".txt" ;
         std::string stringfile101;
         file_name101 >> stringfile101;

  ofstream fout7;

fout7.open(stringfile101);


for (int mu =0; mu < num_mu_steps; mu++) {
for (int QQ =0; QQ < num_histogram_bins; QQ++)
{fout7 << initial_position << " " << pow(10, mu)*mu_step << " " << double(QQ)/double(num_histogram_bins)  << " " << hist_of_bootstrapped_mean_homozygosities[mu][QQ] <<  endl;
 } 
// Here we output mean homozygosity as a function of mu and error bars - mean plus or minus standard deviation of the mean.
  }
  fout7.close();
 




}





inline void print_histogram_of_single_trial_homozygosities(double alpha, double initial_position, int num_trials, int num_time_steps, double rho_inverse, double hist_of_single_trial_homozygosities[][num_histogram_bins]  )
{
char OUTPUTFILE3[50];
  sprintf(OUTPUTFILE3, "histogram_of_single_trial_homozygosities");
  std::stringstream file_name100;
         file_name100 <<  OUTPUTFILE3 << "alpha_value_"<< alpha << "distance_value_" << setw(7) << setfill('0') << initial_position  << "rho_inverse_" << rho_inverse << ".txt" ;
         std::string stringfile100;
         file_name100 >> stringfile100;

  ofstream fout6;

fout6.open(stringfile100);


for (int mu =0; mu < num_mu_steps; mu++) {
for (int QQ =0; QQ < num_histogram_bins; QQ++)
{fout6 << initial_position << " " << pow(10, mu)*mu_step << " " << double(QQ)/double(num_histogram_bins)  << " " << hist_of_single_trial_homozygosities[mu][QQ] <<  endl;
 } 
// Here we output mean homozygosity as a function of mu and error bars - mean plus or minus standard deviation of the mean.
  }
  fout6.close();




}



inline void print_CDF_of_coalescence_times(double alpha, double initial_position, int num_trials, int num_time_steps, double rho_inverse, double CDF_of_DCT[], double lower_CI_CDF[], double upper_CI_CDF[])
{
char OUTPUTFILE22[50];
  sprintf(OUTPUTFILE22, "CDF_of_coalescence_times");
  std::stringstream file_name999;
         file_name999 <<  OUTPUTFILE22 << "alpha_value_"<< alpha << "distance_value_" << 
         setw(7) << setfill('0') << initial_position  << "rho_inverse_" << rho_inverse << ".txt" ;
         std::string stringfile999;
         file_name999 >> stringfile999; 
  ofstream fout55;
fout55.open(stringfile999);
for (int T =0; T < num_cdf_steps; T++) {
fout55  << alpha <<  " " <<  rho_inverse  << " " << exp(double(T)/2.0) << " "
 << initial_position << " " << CDF_of_DCT[T] << " " << lower_CI_CDF[T] <<  " " << upper_CI_CDF[T] << endl;
 // Here we output mean homozygosity as a function of mu and include error bars
}
fout55.close();



}



inline void print_Mean_Homozygosity(double alpha, double initial_position, int num_trials, int num_time_steps, double rho_inverse, double  mean_homozygosity[], double lower_CI[], double upper_CI[])
{
char OUTPUTFILE2[50];
  sprintf(OUTPUTFILE2, "mean_homozygosity_");
  std::stringstream file_name99;
         file_name99 <<  OUTPUTFILE2 << "alpha_value_"<< alpha << "distance_value_" << setw(7) << setfill('0') << 
         initial_position  << "rho_inverse_" << rho_inverse << ".txt" ;
         std::string stringfile99;
         file_name99 >> stringfile99; 
  ofstream fout5;
fout5.open(stringfile99);
for (int mu =0; mu < num_mu_steps; mu++) {
fout5  << alpha <<  " " <<  rho_inverse  << " " << pow(10, mu)*mu_step << " " << initial_position << 
" " << mean_homozygosity[mu] << " " << lower_CI[mu] <<  " " << upper_CI[mu] << endl;
 // Here we output mean homozygosity as a function of mu and include error bars
}
fout5.close();


}






