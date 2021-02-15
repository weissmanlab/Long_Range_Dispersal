#!/bin/sh
c_2D=10.0
coal=1.0
num_trials=500000
num_timesteps=1000

for ALPHA in 0.5 1.0 1.25 1.5 2.0 2.05
do


# varying distance from 1.2^15 to 1.2^30
    ./DISCRETE_2D $ALPHA 0 $num_trials $num_timesteps $c_2D $coal
    for distance in 15.4 18.5 22.2 26.2 31.9 38.4 46.0 55.2 66.2 79.5 95.4 114.5 137.4 164.8 197.8 237.4; do
        distance_dummy=distance
        ./DISCRETE_2D $ALPHA $distance $num_trials $num_timesteps $c_2D $coal
    done


done


cat mean_homozygosity_2D_scale_parameter* > MH_2D_dummy.txt
cat header.txt MH_2D_dummy.txt > MH_2D_discrete.txt
rm MH_2D_dummy.txt
rm mean_homozygosity_2D_scale_parameter*
