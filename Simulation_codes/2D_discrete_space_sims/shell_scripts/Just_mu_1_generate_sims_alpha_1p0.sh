#!/bin/sh
c_2D=10.0
coal=1.0
num_trials=100000000
num_timesteps=10

for ALPHA in 1.0
do


# varying distance from 1.2^15 to 1.2^30
    ./DISCRETE_2D $ALPHA 0 $num_trials $num_timesteps $c_2D $coal
    for distance in 15 19 22 26 32 38 46 55 66 80 95 115 137 165 198 237; do
        distance_dummy=distance
        ./DISCRETE_2D $ALPHA $distance $num_trials $num_timesteps $c_2D $coal
    done


done


cat mean_homozygosity_2D_scale_parameter* > MH_2D_dummy.txt
cat header.txt MH_2D_dummy.txt > MH_2D_discrete_alpha_1p0_Just_mu_1.txt
rm MH_2D_dummy.txt
rm mean_homozygosity_2D_scale_parameter*
