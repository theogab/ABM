#!/bin/bash

for r in $(seq 0.05 0.45 0.5); do 
  for h in $(seq 0.1 0.8 0.9); do
  	for s in $(seq 0.01 0.49 0.5); do
  		Rscript R/ABM_tests.R $r $h $s 1000;
  	done
  done
done

Rscript R/plot_simulations_ABM.R