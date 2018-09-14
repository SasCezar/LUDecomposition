#!/bin/bash
for i in {2..10}; do
    ./main.o >> "sequential_$i.csv"
done
