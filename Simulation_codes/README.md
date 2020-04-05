Note that simulations were run on mac and linux machines only.  The recursive directory creation relies on unix commands which likely won't work on Windows.  Look in the directory Code_without_recursive_directories for a version of the code doesn't create or change directories.  This should work on any OS.



- compilation examples:


- g++  Levy_flight_1D_coalescence_part1_change_dir.cpp -o3 -lgsl -o part1
- g++ Levy_flight_1D_coalescence_part2_change_dir.cpp -o3 -o part2


- Execution examples:

- ./part1 2.05 0 250000 1000 250



- ./part2 2.05 0 250000 1000 1


- Executre ./part1 or ./part2 without any additional arguments to get list of required arguments for each part.

