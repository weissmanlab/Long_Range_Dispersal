#!/bin/bash



for i in 1 2 3 4 5 6 7 8
do
        let var=5**i
	./part1 1.22 $var 10000 1000 20
done

for i in 1 2 3 4 5 6 7 8
do
        let var=5**i
	./part2 1.22 $var 10000 1000 20 .01
done
