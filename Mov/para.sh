#!/usr/bin/env bash

B=(0.5 0.6 0.7 0.8 0.9 1.0)
eps=(0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5)

for i in "${B[@]}"
do
    for j in "${eps[@]}"
    do
        ./main.exe $i $j
    done
    # echo 
done