#!/bin/bash

TIMEOUT=300

# fix 50 seeds for the RNG for reproducibility of results
SEEDS="7347 7945 1788 5178 3923 130 1077 1815 7455 801
4916 5959 3741 596 9770 8351 9936 1482 7252 3152
2201 551 4748 6911 4221 6421 485 9791 572 7642
2592 9420 5852 9092 6528 4826 3497 3132 4321 2274
3988 6254 271 8196 9335 1582 9784 7887 4842 1308"

QUBIT_COUNTS="50 100"
DEPTHS="10 20 30 40 50 60 70"
MIN_WEIGHT=2
MAX_WEIGHT=4
NSAMPLES=1

mkdir -p data
cd data

for qubit in $QUBIT_COUNTS; do
    for depth in $DEPTHS; do
        for seed in $SEEDS; do
            file="pauli_gadget_${qubit}_${depth}_${MIN_WEIGHT}_${MAX_WEIGHT}_${NSAMPLES}_${seed}"
            echo $file
            if [ -f $file ]; then
                echo "EXISTS $file"
            else
                timeout $TIMEOUT ../../quizx/target/release/pauli_gadget_stabrank $qubit $depth $MIN_WEIGHT $MAX_WEIGHT $NSAMPLES $seed
                if [ $? == 124 ]; then echo "TIMEOUT $file"; fi
            fi
        done
    done
done
