#include "part2_same_dir.h"
int main(int argc, char* argv[])
{         if(argc != 6) {cout << "Wrong number of arguments.  Arguments are alpha, initial distance, number of trials, total number of time steps, rho_inverse " << endl; return 0;} 
   
   const double alpha = atof(argv[1]);  // controls power law tail of jump kernel
  double initial_position = atof(argv[2]) ;  // initial signed distance between individuals
const int num_trials = atoi(argv[3]);
const int num_time_steps = atoi(argv[4]);
  const double rho_inverse = atof(argv[5]); 
  int seed = time(0);
  std::mt19937 generator(seed);
  std::uniform_real_distribution<double> uniform_dist(0.0, num_trials);
  double mean_homozygosity[num_mu_steps] = {0}; //probability of two individuals (lineages) being identical given initial seperation and mu
  double CDF_of_DCT[num_cdf_steps] = {0};
double mean_homozygosity_INDIVIDUAL_TRIAL[num_mu_steps] = {0}; // conditional expectation E[exp(-2mut)|path] ] = E[Hom|path]   
double CDF_of_DCT_INDIVIDUAL_TRIAL[num_cdf_steps] = {0};
  double *Contribution_from_each_trial = new double[num_trials];double *Contribution_from_each_trialEXPONENT= new double[num_trials];
  double *dist_of_coalescent_times = new double[num_time_steps];double **List_of_single_trial_homozygosities = new double*[num_mu_steps];
double **List_of_single_trial_CDF = new double*[num_cdf_steps];
for(int i =0; i < num_trials; i++){Contribution_from_each_trial[i] = 0; Contribution_from_each_trialEXPONENT[i] = 0;}
for(int i =0; i < num_time_steps; i++){dist_of_coalescent_times[i] = 0;}
for (int mu = 0; mu < num_mu_steps; mu++) {List_of_single_trial_homozygosities[mu] = new double[num_trials];}
for (int mu = 0; mu < num_mu_steps; mu++) {for(int i =0; i < num_trials; i++){List_of_single_trial_homozygosities[mu][i] = 0; }}
for (int T = 0; T < num_cdf_steps; T++) {List_of_single_trial_CDF[T] = new double[num_trials];}
for (int T = 0; T < num_cdf_steps; T++) {for(int i =0; i < num_trials; i++){List_of_single_trial_CDF[T][i] = 0; }}

  
/*******************************/
//Now that we're in the proper directory with the entrance and exit times of interest, 
//we read them in and calculate the distribution of coalescence times and the mean homozygosity
/**************************************/
for(int trial =0; trial < num_trials; trial++)
 {  Contribution_from_each_trial[trial] = 0;
     Contribution_from_each_trialEXPONENT[trial] = 0;   
       for( int mu = 0; mu < num_mu_steps; mu++){ mean_homozygosity_INDIVIDUAL_TRIAL[mu] =  0;}
     for( int T = 0; T < num_cdf_steps; T++){CDF_of_DCT_INDIVIDUAL_TRIAL[T] =  0;} 
         
          char INPUTFILE_entrance_and_exit_times[50];
  sprintf(INPUTFILE_entrance_and_exit_times, "entrance_and_exit_times");
  std::stringstream file_name_entrance_and_exit_times;
         file_name_entrance_and_exit_times <<  INPUTFILE_entrance_and_exit_times  << 
         "alpha" << alpha << "distance" << initial_position <<  "trial" << trial << ".txt" ;
         std::string stringfile_entrance_and_exit_times;
         file_name_entrance_and_exit_times >> stringfile_entrance_and_exit_times; 
    ifstream fin_entrance_and_exit_times;
    /****************************************/
    // Here we're opening the file containing entrance times associated with each individual trial
    fin_entrance_and_exit_times.open(stringfile_entrance_and_exit_times);
 double entrance_time = -1.0 ;  // If file is empty entrance and exit time will be the same and while loops will be ignored - dist of coalescent times will remain zero
  double exit_time= -1.0 ;

////TIME LOOP BELOW

for (int time =0; time < num_time_steps; time++) {
if (timestep*double(time) >= entrance_time && timestep*double(time) < exit_time) //update exponent and add to dist_of_coalescent_times[time] when in coalescence zone
{Contribution_from_each_trial[trial] =  exp(-Contribution_from_each_trialEXPONENT[trial])*(1- exp(-delta_function_height*rho_inverse*timestep));
 dist_of_coalescent_times[time] += Contribution_from_each_trial[trial]/num_trials;
 Contribution_from_each_trialEXPONENT[trial] += delta_function_height*rho_inverse*timestep;
 // for future times we account for all the time we spent sitting in the coalescence zone during this timestep
for( int mu = 0; mu < num_mu_steps; mu++)
     {mean_homozygosity_INDIVIDUAL_TRIAL[mu] +=  Contribution_from_each_trial[trial]*exp(-pow(10, mu)*2*mu_step*time*timestep); // extra factor of 2 in exponent of Laplace transform is standard in pop gen
      List_of_single_trial_homozygosities[mu][trial] = mean_homozygosity_INDIVIDUAL_TRIAL[mu];   }
for( int T = 0; T < num_cdf_steps; T++)
  {if(exp(double(T)/2.0) > time ){
    CDF_of_DCT_INDIVIDUAL_TRIAL[T] += Contribution_from_each_trial[trial]*timestep;
     List_of_single_trial_CDF[T][trial] = CDF_of_DCT_INDIVIDUAL_TRIAL[T];}}
}
if (timestep*double(time) >= exit_time){fin_entrance_and_exit_times >> entrance_time >> exit_time;}
}
/// TIME LOOP ABOVE
  fin_entrance_and_exit_times.close();
 // bin histogram of single trial homozygosities here.
         for (int mu =0; mu < num_mu_steps; mu++){for(int QQ = 0; QQ < num_histogram_bins; QQ++)
             {if( mean_homozygosity_INDIVIDUAL_TRIAL[mu]  >= double(QQ)/double(num_histogram_bins) && mean_homozygosity_INDIVIDUAL_TRIAL[mu] < double(QQ + 1)/double(num_histogram_bins) )
                {hist_of_single_trial_homozygosities[mu][QQ] += 1/double(num_trials);}
              }}         
           for (int T =0; T < num_cdf_steps; T++){for(int QQ = 0; QQ < num_histogram_bins; QQ++)
             {if( CDF_of_DCT_INDIVIDUAL_TRIAL[T]  >= double(QQ)/double(num_histogram_bins) && CDF_of_DCT_INDIVIDUAL_TRIAL[T] < double(QQ + 1)/double(num_histogram_bins) )
                { hist_of_single_trial_cdfs[T][QQ] += 1/double(num_trials);  }           
              }}

}
//Laplace transform dist of coalescent times to get mean homozygosity
for (int mu =0; mu < num_mu_steps; mu++){for (int time =0; time < num_time_steps; time++) {mean_homozygosity[mu] += dist_of_coalescent_times[time]*exp(-pow(10,mu)*2*mu_step*time*timestep); }}
for (int time =0; time < num_time_steps; time++) { for(int T =0; T < num_cdf_steps; T++) {if(exp(double(T)/2.0) > time){CDF_of_DCT[T] += dist_of_coalescent_times[time]*timestep;}}}
// The loops above calculates the mean homozygosity and CDF



/*******************************************/
// The loops below calculates the confidence interval of the mean homozygosity and CDF
double dummy_rand1 = 0;
int dummy_rand1_INDEX = 0;
// generate list of sample means
for( int mu = 0; mu < num_mu_steps; mu++)
{for(int sample = 0; sample < num_samples_bootstrapped_means; sample++)
   {  for(int trial =0; trial < num_trials; trial++)
         {  dummy_rand1 = floor(uniform_dist(generator));    // we generate random numbers according to the single homozygosity histogram via inverse transform sampling
            dummy_rand1_INDEX = int(dummy_rand1);  //This is the index associated to the above domain value for the Quantile function array
            List_of_bootstrapped_mean_homozygosities[mu][sample] += List_of_single_trial_homozygosities[mu][dummy_rand1_INDEX]/double(num_trials);
          }   
         Sorted_List_of_bootstrapped_mean_homozygosities[mu][sample] = List_of_bootstrapped_mean_homozygosities[mu][sample];

   }
}
// SORT BOOTSTRAPPED HOMOZYGOSITIES BUBBLE SORT 
double tmp1 = 0;
double tmp2 = 0;
for( int mu = 0; mu < num_mu_steps; mu++)
{ for(int sample1 = 0; sample1 < num_samples_bootstrapped_means; sample1++)
    {for(int sample2 = sample1; sample2 < num_samples_bootstrapped_means; sample2++)
    {    tmp1 = Sorted_List_of_bootstrapped_mean_homozygosities[mu][sample1]; tmp2 = Sorted_List_of_bootstrapped_mean_homozygosities[mu][sample2];
         if( tmp2 < tmp1){  Sorted_List_of_bootstrapped_mean_homozygosities[mu][sample1] = tmp2; Sorted_List_of_bootstrapped_mean_homozygosities[mu][sample2] = tmp1;}
    }
    }
 }
// USE SORTED LIST TO GET HISTOGRAM OF BOOSTRAPPED MEAN HOMOZYGOSITIES
//
for( int mu = 0; mu < num_mu_steps; mu++)
{for(int sample = 0; sample < num_samples_bootstrapped_means; sample++)
  { for(int QQ = 0; QQ < num_histogram_bins; QQ++)
   {if( List_of_bootstrapped_mean_homozygosities[mu][sample]  >= double(QQ)/double(num_histogram_bins) && List_of_bootstrapped_mean_homozygosities[mu][sample] < double(QQ + 1)/double(num_histogram_bins) )
                {hist_of_bootstrapped_mean_homozygosities[mu][QQ] += 1/double(num_samples_bootstrapped_means); }
   }  
  }
}


// USE HISTOGRAM TO GET CONDIFENCE INTERVALS
int lower_CI_INDEX[num_mu_steps] = {0};
int upper_CI_INDEX[num_mu_steps] = {0};
for( int mu = 0; mu < num_mu_steps; mu++)
{   lower_CI_INDEX[mu] = int(floor(double(num_samples_bootstrapped_means)*double(.16))); 
    upper_CI_INDEX[mu] = int(floor(double(num_samples_bootstrapped_means)*double(.84)));
}
double lower_CI[num_mu_steps];
double upper_CI[num_mu_steps];
for( int mu = 0; mu < num_mu_steps; mu++)
{lower_CI[mu]= Sorted_List_of_bootstrapped_mean_homozygosities[mu][lower_CI_INDEX[mu]];    // confidence interval for mean homozygosity obtained via bootstrapping
upper_CI[mu]= Sorted_List_of_bootstrapped_mean_homozygosities[mu][upper_CI_INDEX[mu]];
}

// REPEAT THE ABOVE STEPS FOR THE CDF OF COALESCENCE TIMES
dummy_rand1 = 0;
dummy_rand1_INDEX = 0;
// generate list of sample means
for( int T = 0; T < num_cdf_steps; T++)
{for(int sample = 0; sample < num_samples_bootstrapped_cdfs; sample++)
   {  for(int trial =0; trial < num_trials; trial++)
         {dummy_rand1 = floor(uniform_dist(generator));    // we generate random numbers according to the single homozygosity histogram via inverse transform sampling
          dummy_rand1_INDEX = int(dummy_rand1);  //This is the index associated to the above domain value for the Quantile function array
              List_of_bootstrapped_cdfs[T][sample] += List_of_single_trial_CDF[T][dummy_rand1_INDEX]/double(num_trials);
           }       
                Sorted_List_of_bootstrapped_cdfs[T][sample] = List_of_bootstrapped_cdfs[T][sample];
   }

}
// now just sort the list of bootstrapped cdfs via bubble sort
tmp1 = 0;
tmp2 = 0;
for( int T = 0; T < num_cdf_steps; T++)
{ for(int sample1 = 0; sample1 < num_samples_bootstrapped_cdfs; sample1++)
    { for(int sample2 = sample1; sample2 < num_samples_bootstrapped_cdfs; sample2++)
      {  tmp1 = Sorted_List_of_bootstrapped_cdfs[T][sample1];
         tmp2 = Sorted_List_of_bootstrapped_cdfs[T][sample2];
        if( tmp2 < tmp1){  Sorted_List_of_bootstrapped_cdfs[T][sample1] = tmp2; Sorted_List_of_bootstrapped_cdfs[T][sample2] = tmp1; }
      }
    }
 }

// use list of sample cdfs to get histogram of sample cdf values
for( int T = 0; T < num_cdf_steps; T++)
{for(int sample = 0; sample < num_samples_bootstrapped_cdfs; sample++)
  {for(int QQ = 0; QQ < num_histogram_bins; QQ++)
     {  if( List_of_bootstrapped_cdfs[T][sample]  >= double(QQ)/double(num_histogram_bins) && List_of_bootstrapped_cdfs[T][sample] < double(QQ + 1)/double(num_histogram_bins) )
                {  hist_of_bootstrapped_cdfs[T][QQ] += 1/double(num_samples_bootstrapped_means);}
     } 
  }
}

int lower_CI_INDEX_CDF[num_cdf_steps] = {0};
int upper_CI_INDEX_CDF[num_cdf_steps] = {0};
for( int T = 0; T < num_cdf_steps; T++)
{   lower_CI_INDEX_CDF[T] = int(floor(double(num_samples_bootstrapped_cdfs)*double(.16)));  
    upper_CI_INDEX_CDF[T] = int(floor(double(num_samples_bootstrapped_cdfs)*double(.84)));
}
double lower_CI_CDF[num_cdf_steps];
double upper_CI_CDF[num_cdf_steps];
for( int T = 0; T < num_cdf_steps; T++)
{lower_CI_CDF[T]= Sorted_List_of_bootstrapped_cdfs[T][lower_CI_INDEX_CDF[T]];    // confidence interval for mean homozygosity obtained via bootstrapping
upper_CI_CDF[T]= Sorted_List_of_bootstrapped_cdfs[T][upper_CI_INDEX_CDF[T]];
}
/*******************************************/
//Now we output our results to files
/*******************************************/
print_Mean_Homozygosity(alpha, initial_position, num_trials, num_time_steps, rho_inverse, mean_homozygosity, lower_CI, upper_CI);
print_CDF_of_coalescence_times(alpha, initial_position, num_trials, num_time_steps, rho_inverse, CDF_of_DCT,  lower_CI_CDF, upper_CI_CDF);
print_histogram_of_single_trial_homozygosities(alpha, initial_position, num_trials, num_time_steps, rho_inverse, hist_of_single_trial_homozygosities);
print_histogram_of_bootstrapped_mean_homozygosities(alpha, initial_position, num_trials, num_time_steps, rho_inverse, hist_of_bootstrapped_mean_homozygosities);
//print_Sorted_List_of_bootstrapped_mean_homozygosities(alpha, initial_position, num_trials, num_time_steps, rho_inverse, Sorted_List_of_bootstrapped_mean_homozygosities);
//print_list_of_single_trial_homozygosities(alpha, initial_position, num_trials, num_time_steps, rho_inverse, List_of_single_trial_homozygosities);

return 0;
}
