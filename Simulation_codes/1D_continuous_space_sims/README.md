Note that simulations were run on mac and linux machines only.  The recursive directory creation in this code relies on unix commands which likely won't work on Windows.  Look in the directory Code_without_recursive_directories for a version of the code doesn't create or change directories.  This should work on any OS.



- compilation examples:


	g++ Levy_flight_1D_coalescence_part1.cpp -o3 -lgsl -o part1
 
	g++ Levy_flight_1D_coalescence_part2.cpp -o3 -o part2


- Execution examples:

	./part1 1.85 10 1000 1000 250


	./part2 1.85 10 1000 1000 250 1


- Execute ./part1 or ./part2 without any additional arguments to get list of required arguments for each part.

