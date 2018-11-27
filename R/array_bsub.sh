#!/bin/bash

for i in $(seq 0 0.1 0.5); do 
  bsub "sh run_ABM.sh $i";
done