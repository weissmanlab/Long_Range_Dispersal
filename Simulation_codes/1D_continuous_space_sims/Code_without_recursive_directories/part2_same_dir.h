#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <cmath>
#include <random>
#include <iostream>
#include <string>
#include <string.h>
#include <fstream>
#include <sstream>
#include <stdlib.h>     /* atof */
#include <type_traits>
#include <cstdlib>
using namespace std;

std::string input_filename;
const double delta_function_width = 1.0; 
const double delta_function_height = 1.0/delta_function_width; 
const int num_mu_steps = 4;  // number of mu increments 
const int num_cdf_steps = 29;  // make this small if you only care about mean homozygosity and isolation by distance
const double timestep = 1.0; //const double timestep = .1; // for deterministic drift term and dist of coalescence.  for finite t_con this must be the same in part1 and part2
const double mu_step = .001;  
const int num_histogram_bins = 2000;
const int num_samples_bootstrapped_cdfs = 10000;
const int num_samples_bootstrapped_means = 10000;
double hist_of_single_trial_homozygosities[num_mu_steps][num_histogram_bins] = {0};
double hist_of_bootstrapped_mean_homozygosities[num_mu_steps][num_histogram_bins] = {0};
double List_of_bootstrapped_mean_homozygosities[num_mu_steps][num_samples_bootstrapped_means] = {0};
double Sorted_List_of_bootstrapped_mean_homozygosities[num_mu_steps][num_samples_bootstrapped_means] = {0};
double hist_of_single_trial_cdfs[num_cdf_steps][num_histogram_bins] = {0};
double hist_of_bootstrapped_cdfs[num_cdf_steps][num_histogram_bins] = {0};
double List_of_bootstrapped_cdfs[num_cdf_steps][num_samples_bootstrapped_cdfs] = {0};
double Sorted_List_of_bootstrapped_cdfs[num_cdf_steps][num_samples_bootstrapped_cdfs] = {0};
double mean_homozygosity[num_mu_steps] = {0}; //probability of two individuals (lineages) being identical given initial seperation and mu
double CDF_of_DCT[num_cdf_steps] = {0};
double lower_CI[num_mu_steps];
double upper_CI[num_mu_steps];
double lower_CI_CDF[num_cdf_steps];
double upper_CI_CDF[num_cdf_steps];
double time_value;
double mu_value;
double T_value;
double prob_of_no_coalescence_yet;
double Contribution_from_each_trial;
double Contribution_from_each_trialEXPONENT;
int hist_index;


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



inline string return_input_filename(double alpha, double scale_parameter, double initial_position, int num_trials, int num_time_steps, int trial  )
{ char INPUTFILE_entrance_and_exit_times[50];
  sprintf(INPUTFILE_entrance_and_exit_times, "entrance_and_exit_times");
  std::stringstream file_name_entrance_and_exit_times;
  file_name_entrance_and_exit_times <<  INPUTFILE_entrance_and_exit_times << "_scale_parameter_" << scale_parameter << "_alpha" << alpha << "distance" << initial_position <<  "trial" << trial << ".txt" ;
  std::string stringfile_entrance_and_exit_times;
  file_name_entrance_and_exit_times >> stringfile_entrance_and_exit_times; 
  return stringfile_entrance_and_exit_times;
}

double* return_1d_pointer(int index1_size )
{    double *POINTER_ARRAY_1D = new double[index1_size];
     for (int i = 0; i < index1_size; i++) {POINTER_ARRAY_1D[i] = 0;}    
     return POINTER_ARRAY_1D;
}



double** return_2d_pointer(int index1_size, int index2_size )
{    double **POINTER_ARRAY_2D = new double*[index1_size];
     for (int i = 0; i < index1_size; i++) {POINTER_ARRAY_2D[i] = new double[index2_size];}
     for(int i = 0; i < index1_size; i++){for(int j = 0; j < index2_size; j++){POINTER_ARRAY_2D[i][j] = 0; }}
     return POINTER_ARRAY_2D;
}



inline void BOOTSTRAP_CI_MH_and_CDF(double alpha, double initial_position, int num_trials, int num_time_steps, double rho_inverse, double **List_of_single_trial_homozygosities, double **List_of_single_trial_CDF )
{
// SORT BOOTSTRAPPED HOMOZYGOSITIES VIA BUBBLE SORT 

  int seed = 123;
  std::mt19937 generator(seed);
  std::uniform_int_distribution<int> uniform_dist_INT(0.0, num_trials -1);



  int dummy_rand1_INDEX;
// generate list of sample means
  for( int mu = 0; mu < num_mu_steps; mu++){for(int sample = 0; sample < num_samples_bootstrapped_means; sample++)
     {  for(int trial =0; trial < num_trials; trial++)
           {  
              dummy_rand1_INDEX = uniform_dist_INT(generator);
              List_of_bootstrapped_mean_homozygosities[mu][sample] += List_of_single_trial_homozygosities[mu][dummy_rand1_INDEX]/double(num_trials);
            }   
           Sorted_List_of_bootstrapped_mean_homozygosities[mu][sample] = List_of_bootstrapped_mean_homozygosities[mu][sample];
    
      }}



// REPEAT THE ABOVE STEPS FOR THE CDF OF COALESCENCE TIMES
  dummy_rand1_INDEX = 0;
  // generate list of sample means
  for( int T = 0; T < num_cdf_steps; T++)
  {for(int sample = 0; sample < num_samples_bootstrapped_cdfs; sample++)
     {  for(int trial =0; trial < num_trials; trial++)
           {    dummy_rand1_INDEX = uniform_dist_INT(generator);
                List_of_bootstrapped_cdfs[T][sample] += List_of_single_trial_CDF[T][dummy_rand1_INDEX]/double(num_trials);
             }       
            Sorted_List_of_bootstrapped_cdfs[T][sample] = List_of_bootstrapped_cdfs[T][sample];
     }

  }


  double tmp1;
  double tmp2;
  for( int mu = 0; mu < num_mu_steps; mu++)
  { for(int sample1 = 0; sample1 < num_samples_bootstrapped_means; sample1++)
    { for(int sample2 = sample1; sample2 < num_samples_bootstrapped_means; sample2++)
      {    
        tmp1 = Sorted_List_of_bootstrapped_mean_homozygosities[mu][sample1]; 
        tmp2 = Sorted_List_of_bootstrapped_mean_homozygosities[mu][sample2];
        if( tmp2 < tmp1){  Sorted_List_of_bootstrapped_mean_homozygosities[mu][sample1] = tmp2; 
                           Sorted_List_of_bootstrapped_mean_homozygosities[mu][sample2] = tmp1;}
       }}
   }
  // USE SORTED LIST TO GET HISTOGRAM OF BOOSTRAPPED MEAN HOMOZYGOSITIES

  for( int mu = 0; mu < num_mu_steps; mu++){for(int sample = 0; sample < num_samples_bootstrapped_means; sample++)
      {    
         hist_index = int(floor(List_of_bootstrapped_mean_homozygosities[mu][sample]*double(num_histogram_bins)*(1.0 - 1.0e-10))); // 1.0 - 1.0e-10 ensures bin is between 0 and num_histogram_bins -1
         hist_of_bootstrapped_mean_homozygosities[mu][hist_index] += 1/double(num_samples_bootstrapped_means); 
      }}



  for( int mu = 0; mu < num_mu_steps; mu++)
  {   int lower_CI_INDEX = int(floor(double(num_samples_bootstrapped_means)*double(.16))); 
      int upper_CI_INDEX = int(floor(double(num_samples_bootstrapped_means)*double(.84)));
      lower_CI[mu]= Sorted_List_of_bootstrapped_mean_homozygosities[mu][lower_CI_INDEX];    // confidence interval for mean homozygosity obtained via bootstrapping
      upper_CI[mu]= Sorted_List_of_bootstrapped_mean_homozygosities[mu][upper_CI_INDEX];

  }


  for( int T = 0; T < num_cdf_steps; T++)
  { for(int sample1 = 0; sample1 < num_samples_bootstrapped_cdfs; sample1++)
      { for(int sample2 = sample1; sample2 < num_samples_bootstrapped_cdfs; sample2++)
        {  tmp1 = Sorted_List_of_bootstrapped_cdfs[T][sample1];
           tmp2 = Sorted_List_of_bootstrapped_cdfs[T][sample2];
           if( tmp2 < tmp1){  Sorted_List_of_bootstrapped_cdfs[T][sample1] = tmp2; 
                              Sorted_List_of_bootstrapped_cdfs[T][sample2] = tmp1; }
        }
      }
   }

  // use list of sample cdfs to get histogram of sample cdf values
  for( int T = 0; T < num_cdf_steps; T++){for(int sample = 0; sample < num_samples_bootstrapped_cdfs; sample++)
    {hist_index = int(floor(List_of_bootstrapped_cdfs[T][sample]*double(num_histogram_bins)*(1.0 - 1.0e-10)));// 1.0 - 1.0e-10 ensures bin is between 0 and num_histogram_bins -1
     hist_of_bootstrapped_cdfs[T][hist_index] += 1/double(num_samples_bootstrapped_means);
    }}



  for( int T = 0; T < num_cdf_steps; T++)
  {   int lower_CI_INDEX_CDF = int(floor(double(num_samples_bootstrapped_cdfs)*double(.16)));  
      int upper_CI_INDEX_CDF = int(floor(double(num_samples_bootstrapped_cdfs)*double(.84)));
      lower_CI_CDF[T]= Sorted_List_of_bootstrapped_cdfs[T][lower_CI_INDEX_CDF];    // confidence interval for mean homozygosity obtained via bootstrapping
      upper_CI_CDF[T]= Sorted_List_of_bootstrapped_cdfs[T][upper_CI_INDEX_CDF];

  }




}


inline void LAPLACE_TRANSFORM_DCT_TO_MH(double *dist_of_coalescent_times, int num_time_steps)

{

for (int mu =0; mu < num_mu_steps; mu++){for (int time =0; time < num_time_steps; time++) 
  {mean_homozygosity[mu] += dist_of_coalescent_times[time]*exp(-pow(10,mu)*2*mu_step*time*timestep); }}


}



inline void INTEGRATE_DCT_TO_CDF(double *dist_of_coalescent_times, int num_time_steps)

{
            
for(int T =0; T < num_cdf_steps; T++) {for (int time =0; time < num_time_steps; time++) 
  { T_value = exp(double(T)/2.0)*timestep;
     time_value = time*timestep;
    if(T_value > time){CDF_of_DCT[T] += dist_of_coalescent_times[time]*timestep;}

  }}



}


inline void Bin_single_trial_MH_and_CDF_into_histogram(int num_trials, int trial, double** MH_pointer, double** CDF_pointer)

{   for (int mu =0; mu < num_mu_steps; mu++){
              hist_index = int(floor(MH_pointer[mu][trial]*double(num_histogram_bins)*(1.0 - 1.0e-10)));
              hist_of_single_trial_homozygosities[mu][hist_index] += 1/double(num_trials);}                  

    for (int T =0; T < num_cdf_steps; T++){
              hist_index = int(floor(CDF_pointer[T][trial]*double(num_histogram_bins)*(1.0 - 1.0e-10)));
              hist_of_single_trial_cdfs[T][hist_index] += 1/double(num_trials);}


}


inline void print_list_of_single_trial_homozygosities(double alpha, double scale_parameter, double initial_position, int num_trials, int num_time_steps, double rho_inverse, double** List_of_single_trial_homozygosities  )
{
char OUTPUTFILE6[50];
  sprintf(OUTPUTFILE6, "list_of_single_trial_homozygosities");
  std::stringstream file_name103;
         file_name103 <<  OUTPUTFILE6   << "_scale_parameter_" << scale_parameter << "_alpha_value_"<< alpha << "distance_value_" << setw(7) << setfill('0') << initial_position  << "rho_inverse_" << rho_inverse << ".txt" ;
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


inline void print_Sorted_List_of_bootstrapped_mean_homozygosities(double alpha, double scale_parameter, double initial_position, int num_trials, int num_time_steps, double rho_inverse, double Sorted_List_of_bootstrapped_mean_homozygosities[][num_samples_bootstrapped_means]  )
{
char OUTPUTFILE5[50];
  sprintf(OUTPUTFILE5, "sorted_list_of_bootstrapped_mean_homozygosities");
  std::stringstream file_name102;
         file_name102 <<  OUTPUTFILE5  << "_scale_parameter_"<< scale_parameter << "_alpha_value_"<< alpha << "distance_value_" << setw(7) << setfill('0') << initial_position  << "rho_inverse_" << rho_inverse << ".txt" ;
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




inline void print_histogram_of_bootstrapped_mean_homozygosities(double alpha, double scale_parameter, double initial_position, int num_trials, int num_time_steps, double rho_inverse, double hist_of_bootstrapped_mean_homozygosities[][num_histogram_bins]  )
{
char OUTPUTFILE4[50];
  sprintf(OUTPUTFILE4, "histogram_of_bootstrapped_mean_homozygosities");
  std::stringstream file_name101;
         file_name101 <<  OUTPUTFILE4 << "_scale_parameter_"<< scale_parameter  << "_alpha_value_"<< alpha << "distance_value_" << setw(7) << setfill('0') << initial_position  << "rho_inverse_" << rho_inverse << ".txt" ;
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


inline void print_histogram_of_single_trial_homozygosities(double alpha, double scale_parameter, double initial_position, int num_trials, int num_time_steps, double rho_inverse, double hist_of_single_trial_homozygosities[][num_histogram_bins]  )
{
char OUTPUTFILE3[50];
  sprintf(OUTPUTFILE3, "histogram_of_single_trial_homozygosities");
  std::stringstream file_name100;
         file_name100 <<  OUTPUTFILE3 << "_scale_parameter_"<< scale_parameter  << "_alpha_value_"<< alpha << "distance_value_" << setw(7) << setfill('0') << initial_position  << "rho_inverse_" << rho_inverse << ".txt" ;
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

inline void print_CDF_of_coalescence_times(double alpha, double scale_parameter, double initial_position, int num_trials, int num_time_steps, double rho_inverse, double CDF_of_DCT[], double lower_CI_CDF[], double upper_CI_CDF[])
{
  char OUTPUTFILE22[50];
  sprintf(OUTPUTFILE22, "CDF_of_coalescence_times");
  std::stringstream file_name999;
  file_name999 <<  OUTPUTFILE22 << "_scale_parameter_"<< scale_parameter << "_alpha_value_"<< alpha << "distance_value_" << 
  setw(7) << setfill('0') << initial_position  << "rho_inverse_" << rho_inverse << ".txt" ;
  std::string stringfile999;
  file_name999 >> stringfile999; 
  ofstream fout55;
  fout55.open(stringfile999);
  for (int T =0; T < num_cdf_steps; T++) {
    fout55 << scale_parameter << " "  << alpha <<  " " <<  rho_inverse  << " " << exp(double(T)/2.0) << " "
    << initial_position << " " << CDF_of_DCT[T] << " " << lower_CI_CDF[T] <<  " " << upper_CI_CDF[T] << endl;
 // Here we output mean homozygosity as a function of mu and include error bars
   }
   fout55.close();

}


inline void print_Mean_Homozygosity(double alpha, double scale_parameter, double initial_position, int num_trials, int num_time_steps, double rho_inverse, double  mean_homozygosity[], double lower_CI[], double upper_CI[])
{
char OUTPUTFILE2[50];
  sprintf(OUTPUTFILE2, "mean_homozygosity_");
  std::stringstream file_name99;
         file_name99 <<  OUTPUTFILE2  << "scale_parameter_"<< scale_parameter << "_alpha_value_"<< alpha << "distance_value_" << setw(7) << setfill('0') << 
         initial_position  << "rho_inverse_" << rho_inverse << ".txt" ;
         std::string stringfile99;
         file_name99 >> stringfile99; 
  ofstream fout5;
fout5.open(stringfile99);
for (int mu =0; mu < num_mu_steps; mu++) {
fout5 << scale_parameter << " " << alpha <<  " " <<  rho_inverse  << " " << pow(10, mu)*mu_step << " " << initial_position << 
" " << mean_homozygosity[mu] << " " << lower_CI[mu] <<  " " << upper_CI[mu] << endl;
 // Here we output mean homozygosity as a function of mu and include error bars
}
fout5.close();


}
