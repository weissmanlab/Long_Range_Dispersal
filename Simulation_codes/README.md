
Note that simulations were run on mac and linux machines only.  The recursive directory creation relies on unix commands which likely won't work on Windows.  Look in the directory Code_without_recursive_directories for a version of the code doesn't create or change directories.  This should work on any OS.

compilation example:


g++  Levy_flight_1D_coalescence_part1_change_dir.cpp -o3 -lgsl -o part1
g++ Levy_flight_1D_coalescence_part2_change_dir.cpp -o3 -o part2


Execution examples:

./part1 2.05 0 250000 1000 250
./part1 2.05 1 250000 1000 250
./part1 2.05 2.7 250000 1000 250
./part1 2.05 7.39 250000 1000 250
./part1 2.05 20.09 250000 1000 250
./part1 2.05 54.60 250000 1000 250
./part1 2.05 148.41 250000 1000 250
./part1 2.05 403.43 250000 1000 250
./part1 2.05 1096.63 250000 1000 250
./part1 2.05 2980.96 250000 1000 250
./part1 2.05 8103.08 250000 1000 250

./part1 2.05 22026 250000 1000 250

./part1 2.05 59874 1250000 1000 250

./part1 2.05 162755 1250000 1000 250


./part2 2.05 0 250000 1000 1
./part2 2.05 1 250000 1000 1
./part2 2.05 2.7 250000 1000 1
./part2 2.05 7.39 250000 1000 1
./part2 2.05 20.09 250000 1000 1
./part2 2.05 54.60 250000 1000 1
./part2 2.05 148.41 250000 1000 1
./part2 2.05 403.43 250000 1000 1
./part2 2.05 1096.63 250000 1000 1
./part2 2.05 2980.96 250000 1000 1
./part2 2.05 8103.08 250000 1000 1

./part2 2.05 22026 250000 1000 1

./part2 2.05 59874 1250000 1000 1
./part2 2.05 162755 1250000 1000 1




