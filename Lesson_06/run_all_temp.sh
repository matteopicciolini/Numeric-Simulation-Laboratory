#!/bin/bash

L=0.5,2.0  # intervallo da 0.5 a 2.0
N=50  # diviso in 5 parti uguali

start=$(echo $L | awk -F ',' '{print $1}')
end=$(echo $L | awk -F ',' '{print $2}')
step=$(echo "scale=2; ($end - $start) / $N" | bc)


for ((i=0; i<N; i++)); do
    t=$(echo "scale=2; $start + ($i * $step) + ($step / 2)" | bc)
    sed -i '1s/.*/'${t}'/' input-output/input.dat
    ./Exercise_06.1
done