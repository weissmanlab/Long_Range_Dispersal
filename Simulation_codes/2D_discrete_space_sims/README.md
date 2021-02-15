Note that simulations were run on mac and linux machines only. The code should work with any OS and compiler provided that you have the GSL library installed.



- compilation examples:

	g++ DISCRETE_Levy_flight_2D_coalescence.cpp -o3 -lgsl -o DISCRETE_2D

- Execution examples:

	./DISCRETE_2D 1.5 5 1000 1000 10 1

- Execute ./DISCRETE_2D without any additional arguments to get list of required arguments.  

- Parameter values used in simulations can be found in the shell scripts within this directory.

