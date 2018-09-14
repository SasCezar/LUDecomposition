#!/bin/bash
for i in {3..10}; do
    ./main.o >> "omp_$i.csv"
done
