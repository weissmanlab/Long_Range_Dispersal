#!/bin/bash



for i in 1 
do
        let var=5**i
	./part1 .75 0 100000 150 1
done

for i in 1 
do
        let var=5**i
	./part2 .75 0 100000 150 1 .001
done
