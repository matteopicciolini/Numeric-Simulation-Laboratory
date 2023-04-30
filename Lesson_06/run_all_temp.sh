#!/bin/bash

L=0.5,2.0  # intervallo da 0.5 a 2.0
N=50  # numero di punti

start=$(echo $L | awk -F ',' '{print $1}')
end=$(echo $L | awk -F ',' '{print $2}')
step=$(echo "scale=10; ($end - $start) / ($N - 1)" | bc)

for ((i=0; i<N; i++)); do
    t=$(echo "scale=10; $start + ($i * $step)" | bc)
    sed -i "1s/.*/$t/" input-output/input.dat
    echo -e "\n" | ./Exercise_06.1
done