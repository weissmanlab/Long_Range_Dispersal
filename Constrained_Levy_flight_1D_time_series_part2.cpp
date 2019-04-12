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
#include <algorithm>
#include <vector>
using namespace std;
int main(int argc, char* argv[])
{         if(argc != 7) {cout << "Wrong number of arguments.  Arguments are alpha, initial distance, number of trials, total number of time steps,  mutation rate, rho_inverse " << endl; return 0;} 
  
//We declare and initialize relevant variables in this section

  const int num_time_steps = atoi(argv[4]);
  const int num_mu_steps = 1;  // number of mu increments in Laplace time 
  const int num_trials = atoi(argv[3]);
  std::mt19937 generator(time(0));
  std::uniform_real_distribution<double> uniform_dist(0.0, num_trials);
  std::uniform_real_distribution<double> uniform_zero_to_one(0.0, num_trials);
  const int num_distance_steps = 1;   //vary initial seperation exponentially for log plot of mean homozygosity as function of x for fixed mu
  const double periodic_boundary = 10000000; //position constrained between -pb and +pb
  const double timestep = 1.0; //const double timestep = .1; // for deterministic drift term and dist of coalescence.  for finite t_con this must be the same in part1 and part2
  const double MU = atof(argv[5]);
  const double mu_step = MU;  // .0001; //change back later 
  //const double t_con_inverse = .000;//.005; //.5 // (1/tcon) also for determinic drift term
  const double rho_inverse = atof(argv[6]); // .1 ; // (1/rho) is for calculation of expectation over paths. Rho is population density.
  const double delta_function_width = 1.0; //atof(argv[5]);  //this must be same as in part 1
  const double delta_function_height = 1.0/delta_function_width;
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
const int num_histogram_bins = 1e4;
double probability_of_hitting_origin; // this is the probability of obtaining a nonzero homozygosity value for an arbitrary path.  
//If the path hits the coalescence zone at any time it's homozygosity is nonzero.





double **List_of_single_trial_homozygosities = new double*[num_mu_steps];
double **List_of_single_trial_WEIGHTS = new double*[num_mu_steps];

//double **List_of_single_trial_WEIGHTS_CDF = new double*[num_mu_steps]; 
for (int mu = 0; mu < num_mu_steps; mu++) {
  List_of_single_trial_homozygosities[mu] = new double[num_trials];
  List_of_single_trial_WEIGHTS[mu] = new double[num_trials];
  //List_of_single_trial_WEIGHTS_CDF[mu] = new double[num_trials];
}

for (int mu = 0; mu < num_mu_steps; mu++) {
  for(int i =0; i < num_trials; i++)
{List_of_single_trial_homozygosities[mu][i] = 0; 
 List_of_single_trial_WEIGHTS[mu][i] = 0;
 //List_of_single_trial_WEIGHTS_CDF[mu][i] = 0;
}}


double *List_of_rand_uni_dist_nums = new double[num_trials]; // list of random numbers drawn from uniform distribution


double **SORTED_List_of_single_trial_homozygosities = new double*[num_mu_steps];
double **SORTED_List_of_single_trial_WEIGHTS = new double*[num_mu_steps];
double **SORTED_List_of_single_trial_WEIGHTS_CDF = new double*[num_mu_steps];

for (int mu = 0; mu < num_mu_steps; mu++) {
  SORTED_List_of_single_trial_homozygosities[mu] = new double[num_trials];
  SORTED_List_of_single_trial_WEIGHTS[mu] = new double[num_trials];
  SORTED_List_of_single_trial_WEIGHTS_CDF[mu] = new double[num_trials];
}

for (int mu = 0; mu < num_mu_steps; mu++) {
  for(int i =0; i < num_trials; i++)
{SORTED_List_of_single_trial_homozygosities[mu][i] = 0; 
 SORTED_List_of_single_trial_WEIGHTS[mu][i] = 0;
 SORTED_List_of_single_trial_WEIGHTS_CDF[mu][i] = 0;
}}








int single_trial_hist_index;
double hist_of_single_trial_homozygosities[num_mu_steps][num_histogram_bins] = {0};
double CDF_of_single_trial_homozygosities[num_mu_steps][num_histogram_bins] = {0};
double hist_of_bootstrapped_mean_homozygosities[num_mu_steps][num_histogram_bins] = {0};


const int num_samples_bootstrapped_means = 10000;


// YOU MAY NEED POINTER BELOW - check this if errors occur w/ executable
double List_of_bootstrapped_mean_homozygosities[num_mu_steps][num_samples_bootstrapped_means] = {0};
double Sorted_List_of_bootstrapped_mean_homozygosities[num_mu_steps][num_samples_bootstrapped_means] = {0};



for(int i =0; i < num_trials; i++)
{Contribution_from_each_trial[i] = 0; Contribution_from_each_trialEXPONENT[i] = 0;}
 for(int i =0; i < num_time_steps; i++)
{dist_of_coalescent_times[i] = 0;}

  //We declare and initialize relevant variables in the above section


//Next we change to the relevant directory
  

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

//Above we opened/moved into the relevant directory

//Now that we're in the proper directory with the entrance and exit times of interest, 
//we read them in and calculate the distribution of coalescence times and the mean homozygosity



char INPUTFILE_PROB_OF_HITTING_ORIGIN[50]; //we're sampling paths that hit the origin and have nonzero mean homozygosity.  

  sprintf(INPUTFILE_PROB_OF_HITTING_ORIGIN, "probability_of_constrained_trajectory_WEIGHT_");
  
  std::stringstream file_name_PROB_OF_HITTING_ORIGIN;
         file_name_PROB_OF_HITTING_ORIGIN <<  INPUTFILE_PROB_OF_HITTING_ORIGIN  << "alpha" << alpha << "distance" << initial_position  << ".txt" ;
         std::string stringfile_PROB_OF_HITTING_ORIGIN;
         file_name_PROB_OF_HITTING_ORIGIN >> stringfile_PROB_OF_HITTING_ORIGIN; 
    ifstream fin_PROB_OF_HITTING_ORIGIN;

fin_PROB_OF_HITTING_ORIGIN.open(stringfile_PROB_OF_HITTING_ORIGIN);
  
  fin_PROB_OF_HITTING_ORIGIN >>  probability_of_hitting_origin;



fin_PROB_OF_HITTING_ORIGIN.close();


for(int trial =0; trial < num_trials; trial++)
 {  Contribution_from_each_trial[trial] = 0;
     Contribution_from_each_trialEXPONENT[trial] = 0;
          
       for( int mu = 0; mu < num_mu_steps; mu++)
     {
      mean_homozygosity_INDIVIDUAL_TRIAL[mu] =  0;
     
     }


          char INPUTFILE_entrance_and_exit_times[50];

  sprintf(INPUTFILE_entrance_and_exit_times, "entrance_and_exit_times");
  
  std::stringstream file_name_entrance_and_exit_times;
         file_name_entrance_and_exit_times <<  INPUTFILE_entrance_and_exit_times  << "alpha" << alpha << "distance" << initial_position << "MU" << MU <<  "trial" << trial << ".txt" ;
         std::string stringfile_entrance_and_exit_times;
         file_name_entrance_and_exit_times >> stringfile_entrance_and_exit_times; 
    ifstream fin_entrance_and_exit_times;
    
    // Here we're opening the file containing entrance times associated with each individual trial
    char INPUTFILE_trajectory_WEIGHT[50];

  sprintf(INPUTFILE_trajectory_WEIGHT, "trajectory_WEIGHT_");
  
  std::stringstream file_name_trajectory_WEIGHT;
         file_name_trajectory_WEIGHT <<  INPUTFILE_trajectory_WEIGHT  << "alpha" << alpha << "distance" << initial_position <<  "MU" << MU <<  "trial" << trial << ".txt" ;
         std::string stringfile_trajectory_WEIGHT;
         file_name_trajectory_WEIGHT >> stringfile_trajectory_WEIGHT; 
    ifstream fin_trajectory_WEIGHT;
// Here we're opening the file containing the weights associated with each individual trial




    
    fin_trajectory_WEIGHT.open(stringfile_trajectory_WEIGHT);
    fin_entrance_and_exit_times.open(stringfile_entrance_and_exit_times);
 
  double TRAJECTORY_WEIGHT;
 fin_trajectory_WEIGHT >> TRAJECTORY_WEIGHT;
//cout <<  TRAJECTORY_WEIGHT << endl;
 double entrance_time = -1.0 ;  // If file is empty entrance and exit time will be the same and while loops will be ignored - dist of coalescent times will remain zero
  double exit_time= -1.0 ;
  
for (int time =0; time < num_time_steps; time++) {
if (timestep*double(time) >= entrance_time && timestep*double(time) < exit_time) //update exponent and add to dist_of_coalescent_times[time] when in coalescence zone
{Contribution_from_each_trial[trial] =  delta_function_height*rho_inverse*exp(-Contribution_from_each_trialEXPONENT[trial]);
 Contribution_from_each_trialEXPONENT[trial] += delta_function_height*rho_inverse*timestep; //updating for NEXT time step 
 //dist_of_coalescent_times[time] += Contribution_from_each_trial[trial]/num_trials;
 dist_of_coalescent_times[time] += TRAJECTORY_WEIGHT*Contribution_from_each_trial[trial];
for( int mu = 0; mu < num_mu_steps; mu++)
     {mean_homozygosity_INDIVIDUAL_TRIAL[mu] +=  Contribution_from_each_trial[trial]*exp(-pow(10, mu)*2*mu_step*time*timestep); // extra factor of 2 in exponent of Laplace transform is standard in pop gen
      }


}

if (timestep*double(time) >= exit_time)
{fin_entrance_and_exit_times >> entrance_time >> exit_time;}







}
  fin_entrance_and_exit_times.close();
  fin_trajectory_WEIGHT.close();



 // bin histogram of single trial homozygosities here.

         for (int mu =0; mu < num_mu_steps; mu++){

         

            single_trial_hist_index = int(floor(double(num_histogram_bins)*mean_homozygosity_INDIVIDUAL_TRIAL[mu])); // bin values to histogram.
            hist_of_single_trial_homozygosities[mu][single_trial_hist_index] += TRAJECTORY_WEIGHT;

            List_of_single_trial_homozygosities[mu][trial] = mean_homozygosity_INDIVIDUAL_TRIAL[mu];
            List_of_single_trial_WEIGHTS[mu][trial] = TRAJECTORY_WEIGHT;
            mean_homozygosity[mu] += TRAJECTORY_WEIGHT*mean_homozygosity_INDIVIDUAL_TRIAL[mu];
            
            }






}





// The loops above calculates the mean homozygosity


// The loops below calculates the confidence interval of the mean homozygosity




//sort list of individual trial homozygosities and the associated weights

for (int mu =0; mu < num_mu_steps; mu++){

         pair<double, double> *pair_list_of_single_trial_homozygosities_and_weights = new pair<double,double>[num_trials];           

          
          for(int trial = 0; trial < num_trials; trial++)
             
             {  pair_list_of_single_trial_homozygosities_and_weights[trial].first = List_of_single_trial_homozygosities[mu][trial];
                pair_list_of_single_trial_homozygosities_and_weights[trial].second = List_of_single_trial_WEIGHTS[mu][trial];
              //cout << List_of_single_trial_WEIGHTS[mu][trial] << endl;
              }
            
          
             sort(pair_list_of_single_trial_homozygosities_and_weights, pair_list_of_single_trial_homozygosities_and_weights + num_trials);


            for(int trial = 0; trial < num_trials; trial++)
             
             {  //cout << pair_list_of_single_trial_homozygosities_and_weights[trial].first << " " <<  pair_list_of_single_trial_homozygosities_and_weights[trial].second << endl;
              
              SORTED_List_of_single_trial_homozygosities[mu][trial] = pair_list_of_single_trial_homozygosities_and_weights[trial].first;
              SORTED_List_of_single_trial_WEIGHTS[mu][trial] = pair_list_of_single_trial_homozygosities_and_weights[trial].second;
              }



            }






//generate CDF weights of sorted single trial homozygosities


for (int mu =0; mu < num_mu_steps; mu++){

         SORTED_List_of_single_trial_WEIGHTS_CDF[mu][0] =  SORTED_List_of_single_trial_WEIGHTS[mu][0];
          for(int trial = 1; trial < num_trials; trial++)
             
             {   SORTED_List_of_single_trial_WEIGHTS_CDF[mu][trial] = SORTED_List_of_single_trial_WEIGHTS[mu][trial] +   SORTED_List_of_single_trial_WEIGHTS_CDF[mu][trial- 1];
                     
 
              }
            
          
            }


// generate quantile function (inverse of CDF) need to round probs in CDF to fit bins of quantile function array

//int num_quantile_bins = 1e6;

/*
for (int mu =0; mu < num_mu_steps; mu++){
int dummy_cdf_index = 0;
for (int prob_bin = 0; prob_bin < num_quantile_bins; prob_bin++ ){//DEFINE num_quantile_bins
    
  //while(prob_bin > CDF_of_single_trial_homozygosities[mu][dummy_cdf_index] && dummy_cdf_index < num_histogram_bins){dummy_cdf_index++}
   

    //if(prob_bin <= CDF_of_single_trial_homozygosities[mu][dummy_cdf_index] )
  
//CDF_of_single_trial_homozygosities[mu][QQ]
  }


}
*/

// Sort list of single trial homozygosities

 

//const int bootstrapped_mean = 0;




double dummy_rand1 = 0;
int dummy_rand1_INDEX = 0;
//int dummy_prob_hit;
//bool hit;
// generate list of sample means

int running_index = 0;
for(int sample = 0; sample < num_samples_bootstrapped_means; sample++)
   { 


       for(int trial =0; trial < num_trials; trial++)
         { List_of_rand_uni_dist_nums[trial] = uniform_zero_to_one(generator);

         
         }
         sort(List_of_rand_uni_dist_nums, List_of_rand_uni_dist_nums + num_trials);  // list of uniformly distributed random numbers is now sorted in ascending order



      



    for( int mu = 0; mu < num_mu_steps; mu++){
       
       running_index = 0;
      

for(int trial =0; trial < num_trials; trial++)
         { 
               
               if(List_of_rand_uni_dist_nums[trial]  <= SORTED_List_of_single_trial_WEIGHTS_CDF[mu][num_trials])
              {
                while(SORTED_List_of_single_trial_WEIGHTS_CDF[mu][running_index] < List_of_rand_uni_dist_nums[trial] ){running_index++;}
              
         
              // now use the running index to add single trial to sample boostrap mean
              List_of_bootstrapped_mean_homozygosities[mu][sample] += SORTED_List_of_single_trial_homozygosities[mu][running_index]/double(num_trials);
               }
                

                  
                   // Don't add zeros to histogram and bootstrap.  
               //Instead sample the conditonal dist od paths that hit origin and multiply the calculated mean and CI's by probability of hitting orign.  
               //This is a less noisy way of determinig the mean, will have tighter confidence intervals.
                  /*
                else{ List_of_bootstrapped_mean_homozygosities[mu][sample] += 0;
                      // CDF NOT normalized to one.  CDF is normalized to probability of trajectory hitting the origin.
                      // Values above this represent probability that trajectories that never hit the origin and have single trial homozygosity value of zero

                    }
                   */


         }

    /*
    for(int trial =0; trial < num_trials; trial++)
         {   
           //  hit = false;
          //while(hit == false)
          //{
          dummy_rand1 = floor(uniform_dist(generator));    // we generate random numbers according to the single homozygosity histogram via inverse transform sampling
          //cout << dummy_rand1;
          dummy_rand1_INDEX = int(dummy_rand1);  //This is the index associated to the above domain value for the Quantile function array
          
          //dummy_prob_hit = uniform_zero_to_one(generator);
         //  if(dummy_prob_hit < List_of_single_trial_WEIGHTS[mu][trial] ) {hit = true;}   

          //}
          //cout << dummy_rand1_INDEX << endl;
         
              

              List_of_bootstrapped_mean_homozygosities[mu][sample] += List_of_single_trial_homozygosities[mu][dummy_rand1_INDEX]/double(num_trials);
           }   
   */
              
                Sorted_List_of_bootstrapped_mean_homozygosities[mu][sample] = List_of_bootstrapped_mean_homozygosities[mu][sample];

              //cout << List_of_bootstrapped_mean_homozygosities[mu][sample] << endl;
   }


}

// now just sort the list of bootstrapped mean homozygosities via bubble sort (list isn't that big so bubble sort is fine)



double tmp1 = 0;
double tmp2 = 0;

for( int mu = 0; mu < num_mu_steps; mu++)
{ for(int sample1 = 0; sample1 < num_samples_bootstrapped_means; sample1++)
    { 
        for(int sample2 = sample1; sample2 < num_samples_bootstrapped_means; sample2++)
    {    
         tmp1 = Sorted_List_of_bootstrapped_mean_homozygosities[mu][sample1];

         tmp2 = Sorted_List_of_bootstrapped_mean_homozygosities[mu][sample2];
                  if( tmp2 < tmp1)
             {     Sorted_List_of_bootstrapped_mean_homozygosities[mu][sample1] = tmp2;

                   Sorted_List_of_bootstrapped_mean_homozygosities[mu][sample2] = tmp1;

             }

    


    }




    }


 }











// use list of sample means to get histogram of sample mean values
//
for( int mu = 0; mu < num_mu_steps; mu++)
{for(int sample = 0; sample < num_samples_bootstrapped_means; sample++)
  {
      for(int QQ = 0; QQ < num_histogram_bins; QQ++)
   {
          if( List_of_bootstrapped_mean_homozygosities[mu][sample]  >= double(QQ)/double(num_histogram_bins) && List_of_bootstrapped_mean_homozygosities[mu][sample] < double(QQ + 1)/double(num_histogram_bins) )
                {

                  hist_of_bootstrapped_mean_homozygosities[mu][QQ] += 1/double(num_samples_bootstrapped_means);
                }




   }

      
         
     

  }

}



int lower_CI_INDEX[4] = {0};

int upper_CI_INDEX[4] = {0};


for( int mu = 0; mu < num_mu_steps; mu++)
{  
    lower_CI_INDEX[mu] = int(floor(double(num_samples_bootstrapped_means)*double(.16)));
    
    upper_CI_INDEX[mu] = int(floor(double(num_samples_bootstrapped_means)*double(.84)));

}


double lower_CI[4];
double upper_CI[4];


for( int mu = 0; mu < num_mu_steps; mu++)
{ 
lower_CI[mu]= Sorted_List_of_bootstrapped_mean_homozygosities[mu][lower_CI_INDEX[mu]];    // confidence interval for mean homozygosity obtained via bootstrapping
upper_CI[mu]= Sorted_List_of_bootstrapped_mean_homozygosities[mu][upper_CI_INDEX[mu]];
}



// now account for probability of hitting origin.  We have condtional mean and CI's (conditioned on path hitting origin). 
// To get marginal dist of all paths we mutiply by probability of hitting origin.  Simply multiply outputs by probability_of_hitting_origin.



//Now we output our results to files


char OUTPUTFILE[50];
  sprintf(OUTPUTFILE, "dist_of_coalescent_times_");
  std::stringstream file_name;
         file_name <<  OUTPUTFILE << "alpha_value_"<< alpha << "distance_value_" << initial_position <<  "MU" << MU << "rho_inverse_" << rho_inverse << ".txt" ;
         std::string stringfile;
         file_name >> stringfile; 
  ofstream fout4;

fout4.open(stringfile);
for (int time =0; time < num_time_steps; time++) {fout4 << time*timestep << " " << probability_of_hitting_origin*dist_of_coalescent_times[time] << endl;}
fout4.close();


char OUTPUTFILE2[50];
  sprintf(OUTPUTFILE2, "mean_homozygosity_");
  std::stringstream file_name99;
         file_name99 <<  OUTPUTFILE2 << "alpha_value_"<< alpha << "distance_value_" << setw(7) << setfill('0') << initial_position <<  "MU" << MU  << "rho_inverse_" << rho_inverse << ".txt" ;
         std::string stringfile99;
         file_name99 >> stringfile99; 

  ofstream fout5;

double SDOM = 0; // Standard deviation of the mean.  Defined in loop below.

fout5.open(stringfile99);
for (int mu =0; mu < num_mu_steps; mu++) {


fout5 << initial_position << " " << pow(10, mu)*mu_step << " " << probability_of_hitting_origin*mean_homozygosity[mu] << " " << probability_of_hitting_origin*lower_CI[mu] <<  " " << probability_of_hitting_origin*upper_CI[mu] << endl;
 // Here we output mean homozygosity as a function of mu and include error bars
}
fout5.close();


char OUTPUTFILE3[50];
  sprintf(OUTPUTFILE3, "histogram_of_single_trial_homozygosities");
  std::stringstream file_name100;
         file_name100 <<  OUTPUTFILE3 << "alpha_value_"<< alpha << "distance_value_" << setw(7) << setfill('0') << initial_position <<  "MU" << MU  << "rho_inverse_" << rho_inverse << ".txt" ;
         std::string stringfile100;
         file_name100 >> stringfile100;

  ofstream fout6;

fout6.open(stringfile100);


for (int mu =0; mu < num_mu_steps; mu++) {

  fout6 << initial_position << " " << pow(10, mu)*mu_step << " " << double(0)/double(num_histogram_bins)  << " " << probability_of_hitting_origin*hist_of_single_trial_homozygosities[mu][0] + (1-probability_of_hitting_origin) <<  endl;
for (int QQ =1; QQ < num_histogram_bins; QQ++)
{fout6 << initial_position << " " << pow(10, mu)*mu_step << " " << double(QQ)/double(num_histogram_bins)  << " " << probability_of_hitting_origin*hist_of_single_trial_homozygosities[mu][QQ] <<  endl;
 } 
// Here we output mean homozygosity as a function of mu and error bars - mean plus or minus standard deviation of the mean.
  }
  fout6.close();
 



char OUTPUTFILE4[50];
  sprintf(OUTPUTFILE4, "histogram_of_bootstrapped_mean_homozygosities");
  std::stringstream file_name101;
         file_name101 <<  OUTPUTFILE4 << "alpha_value_"<< alpha << "distance_value_" << setw(7) << setfill('0') << initial_position <<  "MU" << MU << "rho_inverse_" << rho_inverse << ".txt" ;
         std::string stringfile101;
         file_name101 >> stringfile101;

  ofstream fout7;

fout7.open(stringfile101);


for (int mu =0; mu < num_mu_steps; mu++) {
fout7 << initial_position << " " << pow(10, mu)*mu_step << " " << double(0)/double(num_histogram_bins)  << " " << probability_of_hitting_origin*hist_of_bootstrapped_mean_homozygosities[mu][0] + (1-probability_of_hitting_origin) <<  endl;
for (int QQ =1; QQ < num_histogram_bins; QQ++)
{fout7 << initial_position << " " << pow(10, mu)*mu_step << " " << double(QQ)/double(num_histogram_bins)  << " " << probability_of_hitting_origin*hist_of_bootstrapped_mean_homozygosities[mu][QQ] <<  endl;
 } 
// Here we output mean homozygosity as a function of mu and error bars - mean plus or minus standard deviation of the mean.
  }
  fout7.close();
 


char OUTPUTFILE5[50];
  sprintf(OUTPUTFILE5, "sorted_list_of_bootstrapped_mean_homozygosities");
  std::stringstream file_name102;
         file_name102 <<  OUTPUTFILE5 << "alpha_value_"<< alpha << "distance_value_" << setw(7) << setfill('0') << initial_position <<  "MU" << MU  << "rho_inverse_" << rho_inverse << ".txt" ;
         std::string stringfile102;
         file_name102 >> stringfile102;

  ofstream fout8;

fout8.open(stringfile102);


for (int mu =0; mu < num_mu_steps; mu++) {
for (int QQ =0; QQ < num_samples_bootstrapped_means; QQ++)
{fout8 << initial_position << " " << pow(10, mu)*mu_step << " " << Sorted_List_of_bootstrapped_mean_homozygosities[mu][QQ]  << " " << probability_of_hitting_origin/double(num_samples_bootstrapped_means)<<  endl;
 } 
// Here we output mean homozygosity as a function of mu and error bars - mean plus or minus standard deviation of the mean.
  }
  fout8.close();



char OUTPUTFILE6[50];
  sprintf(OUTPUTFILE6, "list_of_single_trial_homozygosities");
  std::stringstream file_name103;
         file_name103 <<  OUTPUTFILE6 << "alpha_value_"<< alpha << "distance_value_" << setw(7) << setfill('0') << initial_position  <<  "MU" << MU << "rho_inverse_" << rho_inverse << ".txt" ;
         std::string stringfile103;
         file_name103 >> stringfile103;

  ofstream fout9;

fout9.open(stringfile103);


for (int mu =0; mu < num_mu_steps; mu++) {
for (int QQ =0; QQ < num_trials; QQ++)
{fout9 << initial_position << " " << pow(10, mu)*mu_step << " " << List_of_single_trial_homozygosities[mu][QQ]  << " " << List_of_single_trial_WEIGHTS[mu][QQ]*probability_of_hitting_origin<<  endl;
 } 
// Here we output mean homozygosity as a function of mu and error bars - mean plus or minus standard deviation of the mean.
  }
  fout9.close();






chdir("..");
chdir("..");


return 0;
}





