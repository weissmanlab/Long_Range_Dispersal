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
{  



 char OUTPUTFILE_coal_and_absorption_counts[50];
  sprintf(OUTPUTFILE_coal_and_absorption_counts, "coalescence_and_absorption_counts");
  std::stringstream file_name_coal_and_absorption_counts;
         //file_name_sparse <<  OUTPUTFILE_sparse  << "alpha" << alpha << "distance" << initial_position << "scale_parameter" << scale_parameter << ".txt" ;
         file_name_coal_and_absorption_counts <<  OUTPUTFILE_coal_and_absorption_counts  << ".txt" ;
         std::string stringfile_coal_and_absorption_counts;
         file_name_coal_and_absorption_counts >> stringfile_coal_and_absorption_counts; 
    ofstream fout_coal;
fout_coal.open(stringfile_coal_and_absorption_counts);



 char OUTPUTFILE_coal_points[50];
  sprintf(OUTPUTFILE_coal_points, "coalescence_points");
  std::stringstream file_name_coal_points;
         //file_name_sparse <<  OUTPUTFILE_sparse  << "alpha" << alpha << "distance" << initial_position << "scale_parameter" << scale_parameter << ".txt" ;
         file_name_coal_points <<  OUTPUTFILE_coal_points  << ".txt" ;
         std::string stringfile_coal_points;
         file_name_coal_points >> stringfile_coal_points; 
    ofstream fout_coal_points;
fout_coal_points.open(stringfile_coal_points);
fout_coal_points << "alpha distance scale_parameter trial time pos_1 pos_2" << endl;


 char OUTPUTFILE_sparse[50];
  sprintf(OUTPUTFILE_sparse, "sparse_pair_df");
  std::stringstream file_name_sparse;
         //file_name_sparse <<  OUTPUTFILE_sparse  << "alpha" << alpha << "distance" << initial_position << "scale_parameter" << scale_parameter << ".txt" ;
         file_name_sparse <<  OUTPUTFILE_sparse  << ".txt" ;
         std::string stringfile_sparse;
         file_name_sparse >> stringfile_sparse; 
    ofstream fout_sparse;
   fout_sparse.open(stringfile_sparse);
fout_sparse << "alpha distance scale_parameter trial time pos_1 pos_2" << endl;
int coal_count_array_alpha_p5[3][3];
int absorb_count_array_alpha_p5[3][3];

int coal_count_array_alpha_1p5[3][3];
int absorb_count_array_alpha_1p5[3][3];

int coal_count_array_alpha_2[3][3];
int absorb_count_array_alpha_2[3][3];

for(int i=0; i <3; i++){for(int j= 0; j <3; j++){
 coal_count_array_alpha_p5[i][j] = 0;
absorb_count_array_alpha_p5[i][j] = 0;

coal_count_array_alpha_1p5[i][j] = 0;
absorb_count_array_alpha_1p5[i][j] = 0;

coal_count_array_alpha_2[i][j] = 0;
absorb_count_array_alpha_2[i][j] = 0;


}}


  for( int INITIAL_SEPARATION_LOOP_VAR = 0; INITIAL_SEPARATION_LOOP_VAR < 3; INITIAL_SEPARATION_LOOP_VAR++ ){
 for( int SCALE_PARAMETER_LOOP_VAR = 0; SCALE_PARAMETER_LOOP_VAR < 3; SCALE_PARAMETER_LOOP_VAR++ ){
  for( int ALPHA_LOOP_VAR = 0; ALPHA_LOOP_VAR < 3; ALPHA_LOOP_VAR++ ){
double alpha;
 if(ALPHA_LOOP_VAR == 0){alpha = .5;}
 if(ALPHA_LOOP_VAR == 1){alpha = 1.5;}
  if(ALPHA_LOOP_VAR == 2){alpha = 2.0;}
double scale_parameter;
if(SCALE_PARAMETER_LOOP_VAR == 0){scale_parameter = .0005;}
 if(SCALE_PARAMETER_LOOP_VAR == 1){scale_parameter = .005;}
  if(SCALE_PARAMETER_LOOP_VAR == 2){scale_parameter = 10;}
double initial_position;
  if(INITIAL_SEPARATION_LOOP_VAR == 0){initial_position = 10;}
 if(INITIAL_SEPARATION_LOOP_VAR == 1){initial_position = 1000;}
  if(INITIAL_SEPARATION_LOOP_VAR == 2){initial_position = 10000;}

       //if(argc != 6) {cout << "Wrong number of arguments.  Arguments are alpha, initial distance, number of trials, total number of time steps, and scale parameter." << endl; return 0;}  //, and timescale coarse graining." << endl; return 0;} 
  // include periodic boundaries to ensure dist of coalescent times with jumps converges
 // const int distance_off_set = 0;
  const int num_time_steps = 1000000; //atoi(argv[4]);
 // const int time_scale_coarse_graining = 1; //atoi(argv[5]); //number of steps per (in between) recorded step -this is necessary so that we dont exceed memory requirements with arrays that are too large
  const int num_trials = 50; // atoi(argv[3]);
  //const int num_distance_steps = 1;   //vary initial seperation exponentially for log plot of mean homozygosity as function of x for fixed mu
  const double periodic_boundary = DBL_MAX;//100000000000; //position constrained between -pb and +pb
 
  const double levy_param = 1.00;
  const int SEED = 123;
  double absorbing_wall_position = 999999;
  //const double cutoff = 0; // minimum jump size
  const double timestep = 1.0;//const double timestep = .1; // for deterministic drift term
  const double t_con_inverse = .000;//.005; //.5 // (1/tcon) also for determinic drift term
  const double delta_function_width = 10.0; //atof(argv[5]);;
  //const double alpha = atof(argv[1]);  // controls power law tail of jump kernel
  //const double scale_parameter = 1; //atof(argv[5]);//*pow(timestep, 1.0/alpha); // sets scale of levy alpha stable.  Coalescence zone "delta function" is of width one.  In order to test analytical predictions we want c >> 1.
  // the scale parameter c is related to the generalized diffusion constant D as c =(4D*timestep)^(1/alpha)
  double jump_size_levy = 0;
  double jump_size_fisher = 0;
  double jump_size_levy_RAW = 0;
  double jump_size_fisher_RAW = 0;

  double jump_size_levy_RAW_2 = 0;
  double jump_size_fisher_RAW_2 = 0;
double ORIG_STD = 2*sqrt(alpha*(4*alpha -2)/((2*alpha -2)*(2*alpha -2)*(2*alpha -4))); // This is for F -distribution when alpha > 2
//ORIG_STD is the standard deviation of the fisher distribution.  We rescale so that scale_parameter is the standard deviation.
  const gsl_rng_type * T;
  T = gsl_rng_default;
  gsl_rng* r;
   r = gsl_rng_alloc (T);
    //gsl_rng_set(r, time(0));
   gsl_rng_set(r, SEED);

  //std::default_random_engine generator(time(0));
   //std::mt19937 generator(time(0)); // mersenne twister psuedorandom number generator
   std::mt19937 generator(SEED); // mersenne twister psuedorandom number generator
//std::normal_distribution<double> norm_dist(0.0, 4*D);
 
//std::Levy_alpha_stable_distribution<double> fisher_dist(2*alpha,2*alpha);
double gsl_ran_levy(const gsl_rng *r, double c, double Alpha);

//double dummy;

  //double initial_position = 100 //atof(argv[2]) ;  // initial signed distance between individuals
  double current_position;  // difference between position of lineages!!! 
double current_position_raw_1;
double current_position_raw_2;

//Relevant variables declared and initialized above
/********************************/
 


   current_position = fmod(2*initial_position, periodic_boundary);
   current_position_raw_1 = fmod(-1*initial_position, periodic_boundary);
   current_position_raw_2 = fmod(initial_position, periodic_boundary);



for(int trial =0; trial < num_trials; trial++)
 {   //Contribution_from_each_trial[trial] = 0;
     //C
  bool inside_zone = false; 
   bool inside_zone_new; 
  bool coal_check = false; 
  bool wall_check = false; 

    
    
    for (int time =0; time < num_time_steps; time++) {
     
 double signed_step_size = 0;// sqrt(2*D)*norm_dist(generator) -timestep*t_con_inverse*current_position;
 double signed_step_size_RAW = 0;
 double signed_step_size_RAW_2 = 0;

if(alpha <= 2){ jump_size_levy = gsl_ran_levy(r,  scale_parameter, alpha);}
if(alpha > 2){ jump_size_fisher = gsl_ran_fdist(r, 2*alpha, 2*alpha );}
if(alpha <= 2){ jump_size_levy_RAW = gsl_ran_levy(r,  scale_parameter/2, alpha);}
if(alpha > 2){ jump_size_fisher_RAW = gsl_ran_fdist(r, 2*alpha, 2*alpha );}
if(alpha <= 2){ jump_size_levy_RAW_2 = gsl_ran_levy(r,  scale_parameter/2, alpha);}
if(alpha > 2){ jump_size_fisher_RAW_2 = gsl_ran_fdist(r, 2*alpha, 2*alpha );}


if(alpha <= 2)
{signed_step_size +=  levy_param*jump_size_levy;}
if(alpha > 2){ // F distribution is one sided.  Make it two sided
  if(generator() > generator()) {jump_size_fisher = - jump_size_fisher;} // we want both positive and negative jumps
  signed_step_size += (scale_parameter/ORIG_STD)*jump_size_fisher;}
if(alpha <= 2)
{signed_step_size_RAW +=  levy_param*jump_size_levy_RAW;}
if(alpha > 2){ // F distribution is one sided.  Make it two sided
  if(generator() > generator()) {jump_size_fisher_RAW = - jump_size_fisher_RAW;} // we want both positive and negative jumps
  signed_step_size_RAW += .5*(scale_parameter/ORIG_STD)*jump_size_fisher_RAW;}




  if(alpha <= 2)
{signed_step_size_RAW_2 +=  levy_param*jump_size_levy_RAW_2;}
if(alpha > 2){ // F distribution is one sided.  Make it two sided
  if(generator() > generator()) {jump_size_fisher_RAW_2 = - jump_size_fisher_RAW_2;} // we want both positive and negative jumps
  signed_step_size_RAW_2 += .5*(scale_parameter/ORIG_STD)*jump_size_fisher_RAW_2;}





if(abs(current_position_raw_1) > absorbing_wall_position and  abs(current_position_raw_2) < absorbing_wall_position and wall_check == false and coal_check == false){wall_check = true;

if(current_position_raw_1 < 0){fout_sparse << alpha << " " << initial_position << " " << scale_parameter << " " << trial << " " << time << " " << -1*absorbing_wall_position << " " << current_position_raw_2 << endl;}
if(current_position_raw_1 > 0){fout_sparse << alpha << " " << initial_position << " " << scale_parameter << " " << trial << " " << time << " " << absorbing_wall_position << " " << current_position_raw_2 << endl;}



    if(alpha == .5){absorb_count_array_alpha_p5[SCALE_PARAMETER_LOOP_VAR][INITIAL_SEPARATION_LOOP_VAR] += 1;}
  if(alpha == 1.5){absorb_count_array_alpha_1p5[SCALE_PARAMETER_LOOP_VAR][INITIAL_SEPARATION_LOOP_VAR] += 1;}
  if(alpha == 2) {absorb_count_array_alpha_2[SCALE_PARAMETER_LOOP_VAR][INITIAL_SEPARATION_LOOP_VAR] += 1;}

}



if(abs(current_position_raw_2) > absorbing_wall_position and abs(current_position_raw_1) < absorbing_wall_position and wall_check == false  and coal_check == false){wall_check = true;

 if(current_position_raw_2 < 0){fout_sparse << alpha << " " << initial_position << " " << scale_parameter << " " << trial << " " << time << " " << current_position_raw_1 << " " << -1*absorbing_wall_position << endl;}
if(current_position_raw_2 > 0){fout_sparse << alpha << " " << initial_position << " " << scale_parameter << " " << trial << " " << time << " " << current_position_raw_1 << " " << absorbing_wall_position << endl;}
  
  if(alpha == .5){absorb_count_array_alpha_p5[SCALE_PARAMETER_LOOP_VAR][INITIAL_SEPARATION_LOOP_VAR] += 1;}
  if(alpha == 1.5){absorb_count_array_alpha_1p5[SCALE_PARAMETER_LOOP_VAR][INITIAL_SEPARATION_LOOP_VAR] += 1;}
  if(alpha == 2) {absorb_count_array_alpha_2[SCALE_PARAMETER_LOOP_VAR][INITIAL_SEPARATION_LOOP_VAR] += 1;}




}


if(abs(current_position_raw_2) > absorbing_wall_position and abs(current_position_raw_1) > absorbing_wall_position and wall_check == false  and coal_check == false){wall_check = true;
//fout_sparse << alpha << " " << initial_position << " " << scale_parameter << " " << trial << " " << time << " " << current_position_raw_1 << " " << current_position_raw_1 << endl;


if(current_position_raw_2 < 0 and current_position_raw_1 < 0){fout_sparse << alpha << " " << initial_position << " " << scale_parameter << " " << trial << " " << time << " " << -1*absorbing_wall_position << " " << -1*absorbing_wall_position << endl;}
if(current_position_raw_2 > 0 and current_position_raw_1 < 0){fout_sparse << alpha << " " << initial_position << " " << scale_parameter << " " << trial << " " << time << " " << -1*absorbing_wall_position << " " << absorbing_wall_position << endl;}
if(current_position_raw_2 < 0 and current_position_raw_1 > 0){fout_sparse << alpha << " " << initial_position << " " << scale_parameter << " " << trial << " " << time << " " << absorbing_wall_position << " " << -1*absorbing_wall_position << endl;}
if(current_position_raw_2 > 0 and current_position_raw_1 > 0){fout_sparse << alpha << " " << initial_position << " " << scale_parameter << " " << trial << " " << time << " " << absorbing_wall_position << " " << absorbing_wall_position << endl;}


     if(alpha == .5){absorb_count_array_alpha_p5[SCALE_PARAMETER_LOOP_VAR][INITIAL_SEPARATION_LOOP_VAR] += 1;}
  if(alpha == 1.5){absorb_count_array_alpha_1p5[SCALE_PARAMETER_LOOP_VAR][INITIAL_SEPARATION_LOOP_VAR] += 1;}
  if(alpha == 2) {absorb_count_array_alpha_2[SCALE_PARAMETER_LOOP_VAR][INITIAL_SEPARATION_LOOP_VAR] += 1;}

}





  
 if(abs(current_position) <= delta_function_width/2.0 and wall_check == false and coal_check == false)  // use step function with finite width as replacement for delta function
 { inside_zone_new = true; 
    coal_check = true;
 fout_sparse << alpha << " " << initial_position << " " << scale_parameter << " " << trial << " " << time << " " << current_position_raw_1 << " " << current_position_raw_1 << endl;
 fout_coal_points << alpha << " " << initial_position << " " << scale_parameter << " " << trial << " " << time << " " << current_position_raw_1 << " " << current_position_raw_1 << endl;
  if(alpha == .5){coal_count_array_alpha_p5[SCALE_PARAMETER_LOOP_VAR][INITIAL_SEPARATION_LOOP_VAR] += 1;}
  if(alpha == 1.5){coal_count_array_alpha_1p5[SCALE_PARAMETER_LOOP_VAR][INITIAL_SEPARATION_LOOP_VAR] += 1;}
  if(alpha == 2) {coal_count_array_alpha_2[SCALE_PARAMETER_LOOP_VAR][INITIAL_SEPARATION_LOOP_VAR] += 1;}

 }

 else{ inside_zone_new = false; }




inside_zone = inside_zone_new;

 if(coal_check == false and wall_check == false)
{//fout7 << time << " " << current_position_raw_1 << " " << current_position_raw_2 << endl;

   if(time < 100){fout_sparse << alpha << " " << initial_position << " " << scale_parameter<< " " << trial<< " " << time << " " << current_position_raw_1 << " " << current_position_raw_2 << endl;}
if(100 < time && time < 1000 && time %10 == 0){fout_sparse << alpha << " " << initial_position << " " << scale_parameter<< " " << trial << " " << time << " " << current_position_raw_1 << " " << current_position_raw_2 << endl;}
if(1000 < time && time < 10000 && time %100 == 0){fout_sparse << alpha << " " << initial_position << " "  << scale_parameter << " " << trial << " " << time << " " << current_position_raw_1 << " " << current_position_raw_2 << endl;}
if(10000 < time && time < 100000 && time %1000 == 0){fout_sparse << alpha << " " << initial_position << " " << scale_parameter << " " << trial << " " << time << " " << current_position_raw_1 << " " << current_position_raw_2 << endl;}
if(100000 < time  && time %10000 == 0){fout_sparse << alpha << " " << initial_position << " " << scale_parameter << " " << trial << " " << time << " " << current_position_raw_1 << " " << current_position_raw_2 << endl;}



}



if(coal_check == true or wall_check == true)
{

   if(time < 100){fout_sparse << alpha << " " << initial_position << " " << scale_parameter<< " " << trial<< " " << time << " " << "NA" << " " << "NA" << endl;}
if(100 < time && time < 1000 && time %10 == 0){fout_sparse << alpha << " " << initial_position << " " << scale_parameter << " " << trial << " " << time << " " << "NA" << " " << "NA" << endl;}
if(1000 < time && time < 10000 && time %100 == 0){fout_sparse << alpha << " " << initial_position << " " << scale_parameter << " " << trial << " " << time << " " << "NA" << " " << "NA" << endl;}
if(10000 < time && time < 100000 && time %1000 == 0){fout_sparse << alpha << " " << initial_position << " " << scale_parameter << " " << trial << " " << time << " " << "NA" << " " << "NA" << endl;}
if(100000 < time  && time %10000 == 0){fout_sparse << alpha << " " << initial_position << " " << scale_parameter << " " << trial << " " << time << " " << "NA" << " " << "NA" << endl;}


}
 
 
 current_position_raw_1 = fmod((current_position_raw_1 + signed_step_size_RAW),  periodic_boundary) ; 
 current_position_raw_2 = fmod((current_position_raw_2 + signed_step_size_RAW_2),  periodic_boundary) ; 
  current_position = fmod((current_position_raw_2 - current_position_raw_1),  periodic_boundary) ; 

       };
 
 
  current_position = fmod(2*initial_position, periodic_boundary);
  current_position_raw_1 = fmod(-1*initial_position, periodic_boundary);
   current_position_raw_2 = fmod(initial_position, periodic_boundary);


}
 if(alpha == .5){fout_coal << "alpha " << alpha << " scale parameter " << scale_parameter << " distance " << initial_position <<   " total trials: " <<  num_trials  <<  " number that coalesced: " << coal_count_array_alpha_p5[SCALE_PARAMETER_LOOP_VAR][INITIAL_SEPARATION_LOOP_VAR] << " number that hit absorbing boundary " <<  absorb_count_array_alpha_p5[SCALE_PARAMETER_LOOP_VAR][INITIAL_SEPARATION_LOOP_VAR]  << endl;}
 if(alpha == 1.5){fout_coal << "alpha " << alpha << " scale parameter " << scale_parameter << " distance " << initial_position <<   " total trials: " <<  num_trials  <<  " number that coalesced: " << coal_count_array_alpha_1p5[SCALE_PARAMETER_LOOP_VAR][INITIAL_SEPARATION_LOOP_VAR] << " number that hit absorbing boundary " <<  absorb_count_array_alpha_1p5[SCALE_PARAMETER_LOOP_VAR][INITIAL_SEPARATION_LOOP_VAR]  << endl;}
 if(alpha == 2){fout_coal << "alpha " << alpha << " scale parameter " << scale_parameter << " distance " << initial_position <<   " total trials: " <<  num_trials  <<  " number that coalesced: " << coal_count_array_alpha_2[SCALE_PARAMETER_LOOP_VAR][INITIAL_SEPARATION_LOOP_VAR] << " number that hit absorbing boundary " <<  absorb_count_array_alpha_2[SCALE_PARAMETER_LOOP_VAR][INITIAL_SEPARATION_LOOP_VAR]  << endl;}

}}}
fout_sparse.close();
fout_coal.close();
fout_coal_points.close(); 
  return 0;

}
