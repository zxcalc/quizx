#!/bin/bash

TIMEOUT=300

# fix 100 seeds for the RNG for reproducibility of results
SEEDS="7347 7945 1788 5178 3923 130 1077 1815 7455 801
4916 5959 3741 596 9770 8351 9936 1482 7252 3152
2201 551 4748 6911 4221 6421 485 9791 572 7642
2592 9420 5852 9092 6528 4826 3497 3132 4321 2274
3988 6254 271 8196 9335 1582 9784 7887 4842 1308
1162 700 351 1151 5929 8722 1028 3208 8935 7183
5404 3028 954 8436 2637 4684 8459 510 1714 7113
2733 4596 1299 5057 4005 6811 9948 9318 454 9235
4436 4974 5190 6206 3062 9879 2261 9434 199 7870
5596 9656 6169 9108 4921 3260 9426 7866 7712 2243"

QUBIT_COUNTS="50"
DEPTHS="10 20 30 40 50 60 70"
MIN_WEIGHT=2
MAX_WEIGHT=4
NSAMPLES=3

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
