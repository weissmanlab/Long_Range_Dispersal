#!/bin/sh
c_2D=10.0
coal=1.0
num_trials=10000000
num_timesteps=1000

for ALPHA in 1.25
do


# varying distance from 0 to 12
    ./DISCRETE_2D $ALPHA 0 $num_trials $num_timesteps $c_2D $coal
    for distance in 0 1 2 3 5 7 8 10 12; do
        distance_dummy=distance
        ./DISCRETE_2D $ALPHA $distance $num_trials $num_timesteps $c_2D $coal
    done


done


cat mean_homozygosity_2D_scale_parameter* > MH_2D_dummy.txt
cat header.txt MH_2D_dummy.txt > MH_2D_discrete_alpha_1p25.txt
rm MH_2D_dummy.txt
rm mean_homozygosity_2D_scale_parameter*
