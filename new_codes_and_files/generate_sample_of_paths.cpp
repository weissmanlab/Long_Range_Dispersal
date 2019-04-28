// In this code we sample trajectories that begin at a specified distance from the origin and hit the origin at a random time.  
  /* atof */
#include "generate_sample_of_paths.h"

void generate_sample_of_paths(const double ALPHA, const double INIT_DISTANCE, const int NUM_TRIALS, const int NUM_TIME_STEPS, const double SCALE_PARAMETER, const double MUTATION_RATE)
{      using namespace std;
    
    
  const int num_time_steps = NUM_TIME_STEPS -1;  //the total length of the trajectories considered
  const int num_trials = NUM_TRIALS;  //ACTUAL NUMBER OF TRIALS RUN IS SLIGHLY LESS: Total_num_trials - num_time_steps is true number of trials
  const double periodic_boundary = 10000000; //position constrained between -pb and +pb
  const double levy_param = 1.00;
  const double timestep = 1.0;//const double timestep = .1; // for deterministic drift term
  const double t_con_inverse = .000;//.005; //.5 // (1/tcon) also for determinic drift term
  const double delta_function_width = 1.0; //atof(argv[5]);;
  const double alpha = ALPHA;  // controls power law tail of jump kernel
  const double scale_parameter = SCALE_PARAMETER;//*pow(timestep, 1.0/alpha); // sets scale of levy alpha stable.  Coalescence zone "delta function" is of width one.  In order to test analytical predictions we want c >> 1.
  const double MU = MUTATION_RATE;
  // the scale parameter c is related to the generalized diffusion constant D as c =(4D*timestep)^(1/alpha)
  double initial_position = INIT_DISTANCE;  // initial signed distance between individuals
  double current_position = INIT_DISTANCE;


  const gsl_rng_type * T;
  T = gsl_rng_default;
  gsl_rng* r;
   r = gsl_rng_alloc (T);
   gsl_rng_set(r, time(0));
 

   std::mt19937 generator(time(0)); // mersenne twister psuedorandom number generator

std::normal_distribution<double> normal_dist(0.0, 0.1);

//double gsl_ran_levy(const gsl_rng *r, double c, double Alpha);  // randomly generates numbers according to levy stable dist


/*
double* origin_time_weight_array = new double[num_trials];
double* origin_time_array = new double[num_trials];
//std::vector<int> weights;
    
  double origin_time_weight_array_sum = 0;   // origin time can be any time except first time step and last two timesteps
    for(int i=0; i<num_time_steps -2; ++i) {
        
        //if(i ==0){weights.push_back(0); weight_array[i] = 0; }
        // else{weights.push_back(exp(-2*MU*double(i))); weight_array[i] = exp(-2*MU*double(i));}  // the probability of a particular origin time is related to its relevance at the given mutation rate.  

      if(i ==0){ origin_time_weight_array[i] = 0; }
         else{ origin_time_weight_array[i] = exp(-2*MU*double(i));}  // the probability of a particular origin time is related to its relevance at the given mutation rate.  
            origin_time_weight_array_sum += origin_time_weight_array[i];
    }
  unsigned int* origin_time_draws = new unsigned int[num_time_steps];
 //int* SIGNED_origin_time_draws = new  int[num_time_steps];

gsl_ran_multinomial(r, num_time_steps, num_trials, origin_time_weight_array,  origin_time_draws);

 //for(int i=0; i<num_time_steps; ++i) { SIGNED_origin_time_draws[i] = int(origin_time_draws[i]);  }
*/


/*
int dummy_trial_index = 0;
for(int i=0; i<num_time_steps; ++i) {
    //cout <<  SIGNED_origin_time_draws[i] << endl;
    
    
    for(int j = 0; j < origin_time_draws[i]; j++)   // fill up the origin_time_array with all the origin_time_draws
     {       origin_time_array[dummy_trial_index] = i;
                 
           //cout << origin_time_array[dummy_trial_index] << endl;
          dummy_trial_index++;
      }
      


}

*/
//for(int trial = 0; trial < num_trials; trial++)   
//     {  cout << origin_time_array[trial] << endl;}
//cout << dummy_trial_index << endl;


// origin_time is a time between 1 and (num_timesteps -1).  path can't hit origin at time 0 since it is fixed  by initial condition.
   // (num_timestes -1) is last avalilable index of array of length num_timesteps


  
  //std::discrete_distribution<> draw_origin_time(weights.begin(), weights.end());
  


//delete[] origin_time_weight_array; 



//for(int i =0; i < num_time_steps; i++)
//{ cout << origin_time_draws[i] << endl;}



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
    fout_parameter_file << "number of time steps considered " << num_time_steps << endl;
    fout_parameter_file << "number of trials " << num_trials << endl;
    //fout_parameter_file << "t_con_inverse " << t_con_inverse << endl;
    fout_parameter_file << "scale parameter " << scale_parameter << endl;
    fout_parameter_file << "mutation rate " << MU << endl;
    fout_parameter_file.close();


  

for(int trial =0; trial < num_trials; trial++)
 {  

   char OUTPUTFILE3[50];
  sprintf(OUTPUTFILE3, "entrance_and_exit_times");
  std::stringstream file_name3;
         //file_name3 <<  OUTPUTFILE3  << "alpha" << alpha << "distance" << initial_position << "num_time_steps" << num_time_steps <<  "trial" << total_trial_count << ".txt" ;
         file_name3 <<  OUTPUTFILE3  << "alpha" << alpha << "distance" << initial_position << "MU" << MU <<  "trial" << trial << ".txt" ;
         std::string stringfile3;
         file_name3 >> stringfile3; 
    ofstream fout8;
    fout8.open(stringfile3);
 
  

 char OUTPUTFILE4[50];
  sprintf(OUTPUTFILE4, "Extra_step_jump_size");
  std::stringstream file_name4;
         //file_name4 <<  OUTPUTFILE4  << "alpha" << alpha << "distance" << initial_position << "num_time_steps" << num_time_steps <<  "trial" << total_trial_count << ".txt" ;
         file_name4 <<  OUTPUTFILE4 << "alpha" << alpha << "distance" << initial_position <<  "MU" << MU << "trial" << trial << ".txt" ;


         std::string stringfile4;
         file_name4 >> stringfile4; 
    ofstream fout9;
    fout9.open(stringfile4);
  
  bool inside_zone = false; 
   bool inside_zone_new; 

//total_trial_count += 1;


    double free_trajectory[num_time_steps];  // dummy free trajectory generated will be used to construct constrained trajectory
    //double constrained_trajectory[num_time_steps];  // the constrained trajectory we're interested in generating
    double free_step_list[num_time_steps];  // list of jumps taken at each timestep for free trajectory
    double constrained_step_list[NUM_TIME_STEPS]; // THIS LIST HAS ONE EXTRA STEP list of jumps taken at each timestep for constrained trajectory
 //std::cout << gsl_ran_levy(r,  scale_parameter, alpha) << endl;
    for (int time =0; time < num_time_steps; time++) {
     
 

double signed_step_size = gsl_ran_levy(r,  scale_parameter, alpha);
 
 free_step_list[time] = signed_step_size;  // list of jumps taken at each timestep for free trajectory
 free_trajectory[time] = current_position;
 //constrained_step_list[time] = free_step_list[time]; // list of jumps taken at each timestep for constrained trajectory
 current_position = fmod((current_position + signed_step_size),  periodic_boundary) ; 
 
 //fout9 << free_step_list[time] << endl;
       }
  

//std::uniform_real_distribution<double> uniform_dist_extra(1, num_time_steps + .999 );
int SEED_for_origin_time = int(abs(free_step_list[num_time_steps -1])); 
std::mt19937 ORIGIN_TIME_generator(SEED_for_origin_time);
std::exponential_distribution<double> exponential_dist(MU);
int origin_time = int(floor(exponential_dist(ORIGIN_TIME_generator) + .5) + 1);
while( origin_time >= num_time_steps -2){ origin_time = int(floor(exponential_dist(ORIGIN_TIME_generator) + .5)) + 1;}
// keep drawing until our origin time is in the prescribed range.
// This is fine for MU >= .001,  For smaller mu use GSL multinomial to generate the result of a singel trial, and search through the array to obatin the value.  
//For Mu >= 0.001 This should be faster than GSL multinomial since we wont have a resulting array to search.
/*
const gsl_rng_type * T2;
  T2 = gsl_rng_default;
  gsl_rng* r2;
   r2 = gsl_rng_alloc (T2);
   gsl_rng_set(r2, SEED_for_origin_time);
*/



int SEED_for_big_jump_time = int(abs(free_step_list[num_time_steps -2])); 
 std::mt19937 BIG_JUMP_generator(SEED_for_big_jump_time);

 //int origin_time =  origin_time_array[trial]; //draw_origin_time(generator);  // origin_time is a time between 1 and (num_timesteps -1).  path can't hit origin at time 0 since it is fixed  by initial condition.
   //double prob_of_origin_time = exp(-2*MU*double(origin_time))/origin_time_weight_array_sum;
 std::uniform_int_distribution<int> uniform_dist(0.0, origin_time -1 );  // probability of selecting a given chosen time is 1/(origin_time)
   int big_jump_time = (uniform_dist(BIG_JUMP_generator));  // Time at which we shift the trajectory by a large jump to ensure that the endpoint is the origin
  //double prob_of_big_jump_time = 1/double(origin_time);  //we could choose any int from 0 origin_time

  //double prob_of_selecting_origin_time_and_big_jump_time = prob_of_big_jump_time*prob_of_origin_time; //probability of selecting origin time and big jump time with this method

double EXTRA_inserted_step = -free_trajectory[origin_time]; 
double Gaussian_noise = gsl_ran_gaussian(r, .1);
while(abs(Gaussian_noise) > .8){Gaussian_noise = gsl_ran_gaussian(r, .1); }
std::uniform_real_distribution<double> uni_real_dist(EXTRA_inserted_step - .5,EXTRA_inserted_step + .5);
  for (int time =0; time < big_jump_time; time++) {constrained_step_list[time] = free_step_list[time];}
  constrained_step_list[big_jump_time] = uni_real_dist(generator); //draw from narrow uniform dist centered at EXTRA_inserted_step so we can take ratio of densities later                                    // insert this step into the previously generated trajectory
  for (int time =big_jump_time + 1; time < NUM_TIME_STEPS; time++) {constrained_step_list[time] = free_step_list[time -1];}

  //constrained_step_list[big_jump_time] = free_step_list[big_jump_time] - free_trajectory[origin_time]; // correct constrained list so that it will result in trajectory that ends at the origin


 
//fout9  << free_step_list[big_jump_time] <<  endl; 
//fout9 << constrained_step_list[big_jump_time] << endl;


//fout9 << EXTRA_inserted_step << endl; // mean of uniform dist of inserted steps
fout9 << constrained_step_list[big_jump_time] << endl; // our draw
//fout9 << Gaussian_noise << endl;
//if(Gaussian_noise > .5){std::cout << Gaussian_noise  << std::endl;}
//fout9 << chosen_time << endl;
//fout9 << origin_time << endl;
/*
for (int test_time = 0; test_time < num_time_steps; test_time++)
{
    fout9  << free_step_list[test_time] <<  " " << constrained_step_list[test_time] <<  endl; 

  
  //fout9 << chosen_time << endl;
}
 */


  current_position = fmod(initial_position, periodic_boundary); //reset current position to inital position
     
    



     for (int time =0; time < NUM_TIME_STEPS; time++) {
     
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

// The above process will generate trajectories that begin at the specified point and hit the origin at randomly drawn origin_time.  However, the current sampling method is biased, 
//and we must weight/reweight the trajectories to correct for this bias to make it as if we drew paths from the distributuon of all Levy flights.



fout8.close();
fout9.close();
current_position = fmod(initial_position, periodic_boundary); // reset to inital position

}

chdir("../..");

ofstream fout_dummy;
fout_dummy.open("dummy_file.txt");
fout_dummy.close();
//chdir("..");

  //return 0;
  //delete[] weight_array;
//delete[] origin_time_array;
//delete[] origin_time_draws;

}
