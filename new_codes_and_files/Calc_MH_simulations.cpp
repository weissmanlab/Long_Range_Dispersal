
#include "Calc_MH_simulations.h"

void Calc_MH_simulations(const double ALPHA, const double INIT_DISTANCE, const int NUM_TRIALS, const int NUM_TIME_STEPS, const double MUTATION_RATE, const double RHO_INVERSE)
{        using namespace std;
//We declare and initialize relevant variables in this section

  const gsl_rng_type * T;
  T = gsl_rng_default;      // Will be used in boostrap procedure for CI's
  gsl_rng* r;
   r = gsl_rng_alloc (T);
   gsl_rng_set(r, time(0));




  const int num_time_steps = NUM_TIME_STEPS;
  const int num_mu_steps = 1;  // number of mu increments in Laplace time 
  const int num_trials = NUM_TRIALS;
  std::mt19937 generator(time(0));
  //std::uniform_real_distribution<double> uniform_dist(0.0, num_trials);
  //std::uniform_real_distribution<double> uniform_zero_to_one(0.0, 1);
  const int num_distance_steps = 1;   //vary initial seperation exponentially for log plot of mean homozygosity as function of x for fixed mu
  const double periodic_boundary = 10000000; //position constrained between -pb and +pb
  const double timestep = 1.0; //const double timestep = .1; // for deterministic drift term and dist of coalescence.  for finite t_con this must be the same in part1 and part2
  const double MU = MUTATION_RATE;
  const double mu_step = MU;  
  const double rho_inverse = RHO_INVERSE; // .1 ; // (1/rho) is for calculation of expectation over paths. Rho is population density.
  const double delta_function_width = 1.0; //atof(argv[5]);  //this must be same as in part 1
  const double delta_function_height = 1.0/delta_function_width;
  const double alpha = ALPHA;  // controls power law tail of jump kernel
  double mean_homozygosity[num_mu_steps] = {0}; //probability of two individuals (lineages) being identical given initial seperation and mu
double mean_homozygosity_INDIVIDUAL_TRIAL[num_mu_steps] = {0}; // conditional expectation E[exp(-2mut)|path] ] = E[Hom|path]   
//double second_moment_of_homozygosity_INDIVIDUAL_TRIAL[num_mu_steps] = {0}; // Second moment (conditional) E[exp(-4*mu*t)|path] ] = E[Hom^2|path]
//double conditional_VARIANCE_of_homozygosity_INDIVIDUAL_TRIAL[num_mu_steps] = {0}; // conditional variance Var(exp(-2*mu*t)|path)  = E[Hom^2|path] - E[Hom|path]^2
//double mean_homozygosity_VARIANCE_of_Conditional_Expectations[num_mu_steps] = {0};  // Var(E[Hom| path])
//double mean_homozygosity_expectation_value_of_Conditional_VARIANCE[num_mu_steps] = {0};  // E[Var(Hom| path)]
//double mean_homozygosity_TOTAL_VARIANCE[num_mu_steps] = {0} ;  // Var(exp(-2*mu*t)) = Var(Hom) = Var(E[Hom| path]) + E[Var(Hom| path)]
double initial_position = INIT_DISTANCE ;  // initial signed distance between individuals
  double *Contribution_from_each_trial = new double[num_trials];// {0};
  double *Contribution_from_each_trialEXPONENT= new double[num_trials];
  double *dist_of_coalescent_times = new double[num_time_steps];
for(int i =0; i < num_trials; i++)
{Contribution_from_each_trial[i] = 0; Contribution_from_each_trialEXPONENT[i] = 0;}
 for(int i =0; i < num_time_steps; i++)
{dist_of_coalescent_times[i] = 0;}

const int num_histogram_bins = 1e4;
//double probability_of_hitting_origin; // this is the probability of obtaining a nonzero homozygosity value for an arbitrary path.  
//If the path hits the coalescence zone at any time it's homozygosity is nonzero.
int single_trial_hist_index;
double hist_of_single_trial_homozygosities[num_mu_steps][num_histogram_bins] = {0};
double CDF_of_single_trial_homozygosities[num_mu_steps][num_histogram_bins] = {0};
double hist_of_bootstrapped_mean_homozygosities[num_mu_steps][num_histogram_bins] = {0};


const int num_samples_bootstrapped_means = 10000;


// YOU MAY NEED POINTER BELOW - check this if errors occur w/ executable
double List_of_bootstrapped_mean_homozygosities[num_mu_steps][num_samples_bootstrapped_means] = {0};
double Sorted_List_of_bootstrapped_mean_homozygosities[num_mu_steps][num_samples_bootstrapped_means] = {0};


double **List_of_single_trial_homozygosities = new double*[num_mu_steps];
double **List_of_single_trial_WEIGHTS = new double*[num_mu_steps];
double *List_of_rand_uni_dist_nums = new double[num_trials]; // list of random numbers drawn from uniform distribution


double **SORTED_List_of_single_trial_homozygosities = new double*[num_mu_steps];
double **SORTED_List_of_single_trial_WEIGHTS = new double*[num_mu_steps];
double **SORTED_List_of_single_trial_WEIGHTS_CDF = new double*[num_mu_steps];
//double **List_of_single_trial_WEIGHTS_CDF = new double*[num_mu_steps]; 
for (int mu = 0; mu < num_mu_steps; mu++) {
  List_of_single_trial_homozygosities[mu] = new double[num_trials];
  List_of_single_trial_WEIGHTS[mu] = new double[num_trials];
 SORTED_List_of_single_trial_homozygosities[mu] = new double[num_trials];
  SORTED_List_of_single_trial_WEIGHTS[mu] = new double[num_trials];
  SORTED_List_of_single_trial_WEIGHTS_CDF[mu] = new double[num_trials];
  for(int i =0; i < num_trials; i++)
  {List_of_single_trial_homozygosities[mu][i] = 0; 
  List_of_single_trial_WEIGHTS[mu][i] = 0;
  SORTED_List_of_single_trial_homozygosities[mu][i] = 0; 
  SORTED_List_of_single_trial_WEIGHTS[mu][i] = 0;
  SORTED_List_of_single_trial_WEIGHTS_CDF[mu][i] = 0;
  }
}



double* trajectory_weight_array = new double[NUM_TRIALS];
double trajectory_weight_sum = 0;
  //We declare and initialize relevant variables in the above section


//Next we change to the relevant directory
  

  std::stringstream file_name2;
         file_name2 <<  "alpha_value_"  <<  alpha  ;   // This is the directory name
         std::string stringfile2;
         file_name2 >> stringfile2; 
         
       

char OUTPUTFILE_Parent_Directory[50];
strcpy(OUTPUTFILE_Parent_Directory, stringfile2.c_str());


chdir(OUTPUTFILE_Parent_Directory);

     initial_position = INIT_DISTANCE; 
   
      std::stringstream file_name4;
         file_name4 <<  "distance_value_"  <<  initial_position  ;   // This is the directory name
         std::string stringfile4;
         file_name4 >> stringfile4; 
         

char OUTPUTFILE_Child_Directory[50];
strcpy(OUTPUTFILE_Child_Directory, stringfile4.c_str());


chdir(OUTPUTFILE_Child_Directory);












/*
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
*/

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
  //cout << TRAJECTORY_WEIGHT << endl;
  trajectory_weight_array[trial] = TRAJECTORY_WEIGHT;
  trajectory_weight_sum += trajectory_weight_array[trial];
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










//gsl_ran_multinomial (const gsl_rng * r, size_t K, unsigned int N, const double p[], unsigned int n[]);

//sort list of individual trial homozygosities and the associated weights


/*

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



*/


//double dummy_rand1 = 0;
//int dummy_rand1_INDEX = 0;
//int dummy_prob_hit;
//bool hit;
// generate list of sample means

//int running_index = 0;

// Here we sample histogram directly/explicitly

/*

std::uniform_real_distribution<double> uniform_zero_to_one(0.0, 1);
double check_if_nonzero = gsl_ran_binomial(r, trajectory_weight_sum, 1);  // draw yes or no from probability of a nonzero trajectory

for( int mu = 0; mu < num_mu_steps; mu++){

for(int sample = 0; sample < num_samples_bootstrapped_means; sample++)
   {   unsigned int* trajectory_draws = new unsigned int[NUM_TRIALS];
        
   
        gsl_ran_multinomial(r, num_trials, num_trials, trajectory_weight_array,  trajectory_draws);
//double sample_trajectory_weight_sum = 0;
         //int dummy_trial_index = 0;
for(int i=0; i<num_trials; ++i) {
   
    for(int j = 0; j < trajectory_draws[i]; j++)   // fill up the origin_time_array with all the origin_time_draws
     {      check_if_nonzero = gsl_ran_binomial(r, trajectory_weight_sum, 1);  
         if(check_if_nonzero == 1)
         {List_of_bootstrapped_mean_homozygosities[mu][sample] += List_of_single_trial_homozygosities[mu][i]/num_trials;}
                   // total proabability of not hitting zero for this bootstrapped sample

                 //rather than add add an entry with value zero and weight (1-trajectory_weight_sum) we can just multiply all our nonzero values by the probability of a nonzero value, trajectory_weight_sum. 
           //cout << origin_time_array[dummy_trial_index] << endl;
          //dummy_trial_index++;
      }
      //cout << trajectory_draws[i] << endl;
    }
      
    
     
      Sorted_List_of_bootstrapped_mean_homozygosities[mu][sample] = List_of_bootstrapped_mean_homozygosities[mu][sample];
      
      //all the undersampling that comes with  explicitly including zero in the draws.  //we avoid this here

}

}
       
*/

//BELOW GIVES THE CORRECT BOOTSTRAP FOR OUR IMPORTANCE SAMPLING PRODECURE.  THIS IS THE ERROR IN OUR ACTUAL ESTIMATES



std::uniform_int_distribution<int> uniform_dist(0.0, num_trials -1 ); 

int dummy_index_bootstrap;
for( int mu = 0; mu < num_mu_steps; mu++){

for(int sample = 0; sample < num_samples_bootstrapped_means; sample++)
   {   //unsigned int* trajectory_draws = new unsigned int[NUM_TRIALS];
         
   
        //gsl_ran_multinomial(r, num_trials, num_trials, trajectory_weight_array,  trajectory_draws);


//int sample_trajectory_weight_sum = 0;
         //int dummy_trial_index = 0;
for(int trial=0; trial<num_trials; trial++) {
   
        
              dummy_index_bootstrap = uniform_dist(generator);
      List_of_bootstrapped_mean_homozygosities[mu][sample] +=  List_of_single_trial_WEIGHTS[mu][dummy_index_bootstrap]*List_of_single_trial_homozygosities[mu][dummy_index_bootstrap]/num_trials;
                 //sample_trajectory_weight_sum += List_of_single_trial_WEIGHTS[mu][i];  // total proabability of not hitting zero for this bootstrapped sample

                 //rather than add add an entry with value zero and weight (1-trajectory_weight_sum) we can just multiply all our nonzero values by the probability of a nonzero value, trajectory_weight_sum. 
           //cout << origin_time_array[dummy_trial_index] << endl;
          //dummy_trial_index++;
      
      //cout << trajectory_draws[i] << endl;
    }
      
    
     
      Sorted_List_of_bootstrapped_mean_homozygosities[mu][sample] = List_of_bootstrapped_mean_homozygosities[mu][sample];
      // We draw num_trial nonzero draws from histogram.  For our sample we sum the weights of all the draws and multiply by this factor.  
      //This sum is our estimate of the probability of a nonzero draw.  This allows us to estimate the confidence intervals of the MH without dealing with 
      //all the undersampling that comes with  explicitly including zero in the draws

}}









/*



for( int mu = 0; mu < num_mu_steps; mu++){

for(int sample = 0; sample < num_samples_bootstrapped_means; sample++)
   {   unsigned int* trajectory_draws = new unsigned int[NUM_TRIALS];
        
   
        gsl_ran_multinomial(r, num_trials, num_trials, trajectory_weight_array,  trajectory_draws);
//double sample_trajectory_weight_sum = 0;
         //int dummy_trial_index = 0;
for(int i=0; i<num_trials; ++i) {
   
    for(int j = 0; j < trajectory_draws[i]; j++)   // fill up the origin_time_array with all the origin_time_draws
     {     
         List_of_bootstrapped_mean_homozygosities[mu][sample] += trajectory_weight_sum*List_of_single_trial_homozygosities[mu][i]/num_trials;
                   // total proabability of not hitting zero is trajectory weight sum.  This accounts for the fact that the histogram of trials is not normalized
                   // since their is some probability of a trajectory never hitting the origin and producing a single trial homozygosity of zero

                 //rather than add add an entry with value zero and weight (1-trajectory_weight_sum) we can just multiply all our nonzero values by the probability of a nonzero value, trajectory_weight_sum. 
           //cout << origin_time_array[dummy_trial_index] << endl;
          //dummy_trial_index++;
      }
      //cout << trajectory_draws[i] << endl;
    }
      
    
     
      Sorted_List_of_bootstrapped_mean_homozygosities[mu][sample] = List_of_bootstrapped_mean_homozygosities[mu][sample];
      
      //all the undersampling that comes with  explicitly including zero in the draws.  //we avoid this here

}

}



*/



/*


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
               //Instead sample the conditonal dist of paths that hit origin and multiply the calculated mean and CI's by probability of hitting orign.  
               //This is a less noisy way of determinig the confidence intervals.
                  /*
                else{ List_of_bootstrapped_mean_homozygosities[mu][sample] += 0;
                      // CDF NOT normalized to one.  CDF is normalized to probability of trajectory hitting the origin.
                      // Values above this represent probability that trajectories that never hit the origin and have single trial homozygosity value of zero

                    }
                  */


        // }

   
          //      Sorted_List_of_bootstrapped_mean_homozygosities[mu][sample] = List_of_bootstrapped_mean_homozygosities[mu][sample];

              //cout << List_of_bootstrapped_mean_homozygosities[mu][sample] << endl;
//   }


//}




// now just sort the list of bootstrapped mean homozygosities via bubble sort (list isn't that big so bubble sort is fine)



//double tmp1 = 0;
//double tmp2 = 0;

for( int mu = 0; mu < num_mu_steps; mu++){ 


  //for(int sample1 = 0; sample1 < num_samples_bootstrapped_means; sample1++)
    //{ 

       sort(Sorted_List_of_bootstrapped_mean_homozygosities[mu], Sorted_List_of_bootstrapped_mean_homozygosities[mu] + num_samples_bootstrapped_means);
     /*

      for(int sample2 = sample1; sample2 < num_samples_bootstrapped_means; sample2++)
      {    
         tmp1 = Sorted_List_of_bootstrapped_mean_homozygosities[mu][sample1];

         tmp2 = Sorted_List_of_bootstrapped_mean_homozygosities[mu][sample2];
                  if( tmp2 < tmp1)
             {     Sorted_List_of_bootstrapped_mean_homozygosities[mu][sample1] = tmp2;

                   Sorted_List_of_bootstrapped_mean_homozygosities[mu][sample2] = tmp1;

               }

    


       }
  
      */

 // }


}


//for(int sample = 0; sample < num_samples_bootstrapped_means; sample++){cout << Sorted_List_of_bootstrapped_mean_homozygosities[0][sample] << endl;}


// use list of sample means to get histogram of sample mean values
//


/*

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

*/

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



// CONVERT FROM CONDITIONAL TO MARGINAL DISTRIBUTION



////DEBUGGING ONLY - COMMENT BELOW  OUT !!!!!!!!

//probability_of_hitting_origin = 1;  // THIS IS FOR DEBUGGING/FREE PATH COMPARSION.  TURN OFF LATER!!!!! Probl of hitting origin is not 1!!!!!!!!!!!!
////DEBUGGING ONLY - COMMENT ABOVE  OUT !!!!!!!!
/*
for( int mu = 0; mu < num_mu_steps; mu++)
{ 
mean_homozygosity[mu] *= probability_of_hitting_origin;
lower_CI[mu] *= probability_of_hitting_origin;    
upper_CI[mu] *= probability_of_hitting_origin;
}

for (int time =0; time < num_time_steps; time++) { dist_of_coalescent_times[time] *= probability_of_hitting_origin;}
*/




//OUTPUT TO FILES


char OUTPUTFILE[50];
  sprintf(OUTPUTFILE, "dist_of_coalescent_times_");
  std::stringstream file_name;
         file_name <<  OUTPUTFILE << "alpha_value_"<< alpha << "distance_value_" << initial_position <<  "MU" << MU << "rho_inverse_" << rho_inverse << ".txt" ;
         std::string stringfile;
         file_name >> stringfile; 
  ofstream fout4;

fout4.open(stringfile);
for (int time =0; time < num_time_steps; time++) {fout4 << time*timestep << " " << dist_of_coalescent_times[time] << endl;}
fout4.close();


char OUTPUTFILE2[50];
  sprintf(OUTPUTFILE2, "SIMULATION_mean_homozygosity_");
  std::stringstream file_name99;
         file_name99 <<  OUTPUTFILE2 << "alpha_value_" << alpha << "distance_value_"  << initial_position <<  "MU" << MU  << "rho_inverse_" << rho_inverse << ".txt" ;
         std::string stringfile99;
         file_name99 >> stringfile99; 

  ofstream fout5;

double SDOM = 0; // Standard deviation of the mean.  Defined in loop below.

fout5.open(stringfile99);
for (int mu =0; mu < num_mu_steps; mu++) {


fout5 << alpha << " " << rho_inverse <<  " " << MU << " " << initial_position << " "  << mean_homozygosity[mu] << " " << lower_CI[mu] <<  " " << upper_CI[mu] << endl;
 // Here we output mean homozygosity as a function of mu and include error bars
}
fout5.close();








//EXTRA FILES GENERATED BELOW
//Output Extra files below to help check for errors and debug

/*
char OUTPUTFILE3[50];
  sprintf(OUTPUTFILE3, "histogram_of_single_trial_homozygosities");
  std::stringstream file_name100;
         file_name100 <<  OUTPUTFILE3 << "alpha_value_"<< alpha << "distance_value_" << setw(7) << setfill('0') << initial_position <<  "MU" << MU  << "rho_inverse_" << rho_inverse << ".txt" ;
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
 

*/
/*
char OUTPUTFILE4[50];
  sprintf(OUTPUTFILE4, "histogram_of_bootstrapped_mean_homozygosities");
  std::stringstream file_name101;
         file_name101 <<  OUTPUTFILE4 << "alpha_value_"<< alpha << "distance_value_" << setw(7) << setfill('0') << initial_position <<  "MU" << MU << "rho_inverse_" << rho_inverse << ".txt" ;
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
 */


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
{fout8 << initial_position << " " << pow(10, mu)*mu_step << " " << Sorted_List_of_bootstrapped_mean_homozygosities[mu][QQ] <<  endl;
 } 
// Here we output mean homozygosity as a function of mu and error bars - mean plus or minus standard deviation of the mean.
  }
  fout8.close();



char OUTPUTFILE6[50];
  sprintf(OUTPUTFILE6, "list_of_single_trial_homozygosities_and_weights");
  std::stringstream file_name103;
         file_name103 <<  OUTPUTFILE6 << "alpha_value_"<< alpha << "distance_value_" << setw(7) << setfill('0') << initial_position  <<  "MU" << MU << "rho_inverse_" << rho_inverse << ".txt" ;
         std::string stringfile103;
         file_name103 >> stringfile103;

  ofstream fout9;

fout9.open(stringfile103);


for (int mu =0; mu < num_mu_steps; mu++) {
for (int QQ =0; QQ < num_trials; QQ++)
{fout9 << initial_position << " " << pow(10, mu)*mu_step << " " << List_of_single_trial_homozygosities[mu][QQ]  << " "  << List_of_single_trial_WEIGHTS[mu][QQ] << endl;  //<< " " << List_of_single_trial_WEIGHTS[mu][QQ]*probability_of_hitting_origin <<  endl;
 } 
// Here we output mean homozygosity as a function of mu and error bars - mean plus or minus standard deviation of the mean.
  }
  fout9.close();


/*
ofstream fout_rand_nums;


fout_rand_nums.open("sorted_list_of_uni_distrandom_numbers.txt");



for (int trial =0; trial < num_trials; trial++)
{fout_rand_nums <<  List_of_rand_uni_dist_nums[trial] << endl;  //<< " " << List_of_single_trial_WEIGHTS[mu][QQ]*probability_of_hitting_origin <<  endl;
 } 
// Here we output mean homozygosity as a function of mu and error bars - mean plus or minus standard deviation of the mean.
  
  fout_rand_nums.close();
*/



chdir("..");

chdir("..");

//Now we create and move into a directory for storing plots
  std::stringstream file_name_chdir_plot;
         file_name_chdir_plot <<  "MH_plots" ;   // This is the directory name
         std::string stringfile_chdir_plot;
         file_name_chdir_plot >> stringfile_chdir_plot; 
         
          std::stringstream file_name_mkdir_plot;
         file_name_mkdir_plot <<  "mkdir MH_plots"  ;   // This is mkdir directory name
         std::string stringfile_mkdir_plot;
         getline(file_name_mkdir_plot, stringfile_mkdir_plot); 

char Create_Plot_Directory[50];
strcpy(Create_Plot_Directory, stringfile_mkdir_plot.c_str());
system(Create_Plot_Directory);

char Change_to_Plot_Directory[50];
strcpy(Change_to_Plot_Directory, stringfile_chdir_plot.c_str());
chdir(Change_to_Plot_Directory);




char OUTPUTFILE_MH_for_plots[100];
  sprintf(OUTPUTFILE_MH_for_plots, "SIMULATION_mean_homozygosity_");
  std::stringstream file_name_for_mh_plots;
         file_name_for_mh_plots <<  OUTPUTFILE_MH_for_plots << "alpha_value_"<< alpha << "distance_value_"  << initial_position <<  "MU" << MU  << "rho_inverse_" << rho_inverse << ".txt" ;
         std::string stringfile_for_mh_plots;
         file_name_for_mh_plots >> stringfile_for_mh_plots; 

  ofstream fout_for_mh_plots;



fout_for_mh_plots.open(stringfile_for_mh_plots);
for (int mu =0; mu < num_mu_steps; mu++) {


//fout_for_mh_plots << initial_position << " " << pow(10, mu)*mu_step << " " << mean_homozygosity[mu] << " " << lower_CI[mu] <<  " " << upper_CI[mu] << endl;
  fout_for_mh_plots << ALPHA << " "  <<  RHO_INVERSE << " " << MU <<  " " << initial_position << " " << mean_homozygosity[mu] << " " << lower_CI[mu] <<  " " << upper_CI[mu] << endl;
//ADD ALPHA, RHO INVERSE and "SIMULATION" DESIGNATION TO OUTPUT SO THAT ALL OF THESE FILES CAN BE CONCATONATED INTO DATAFRAME FOR PLOTTING

 // Here we output mean homozygosity as a function of mu and include error bars
}
fout_for_mh_plots.close();







chdir("..");


}





