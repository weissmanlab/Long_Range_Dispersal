#!/usr/bin/env python3

import numpy as np

entrance_exit_file = 't_vals.txt' # file storing entrance and exit times to coalescence zone. each line should be in form 'entrance_time exit_time'
rho = [10, 20] # densities to use in calculations
mu = np.logspace(-5, -1) # mutation rates to use in calculations

hom = np.zeros((len(rho), len(mu))) #homozygosities; one entry per (rho, mu) pair
coal_time = 0 # total time spent in coalescence zone

with open(entrance_exit_file, 'r') as infile:
	for line in infile:
		t = [float(x) for x in line.split()] # t = [entrance_time, exit_time]
		for i, n in enumerate(rho):
			for j, u in enumerate(mu):
				hom[i, j] += np.exp(-2 * u * t[0] - coal_time / n) * (-np.expm1(-(2 * u + 1 / n) * np.diff(t))) / (1 + 2 * u * n)
		coal_time += np.diff(t)
				
print(hom)

# upper bound on how much error in homozygosity is introduced by cutting off the trajectory:
error_bound = np.empty((len(rho), len(mu)))
for i, n in enumerate(rho):
	for j, u in enumerate(mu):
		error_bound[i,j] = np.exp(-2 * u * t[1] - coal_time / n)
print(error_bound)