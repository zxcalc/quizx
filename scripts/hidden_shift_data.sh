#!/bin/bash

TIMEOUT=600

# fix 100 seeds for the RNG for reproducibility of results
SEEDS="1337
42
106"

qubit=50
ccz=5

mkdir -p data
cd data

for seed in $SEEDS
do
    if [ -f "hidden_shift_${qubit}_${ccz}_${seed}" ]; then
        echo "EXISTS $qubit $ccz $seed"
    else
        timeout $TIMEOUT ../../quizx/target/release/hidden_shift_stabrank $qubit $ccz $seed
        if [ $? == 124 ]; then echo "TIMEOUT $qubit $ccz $seed"; fi
    fi
done
