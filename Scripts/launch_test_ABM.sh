#!/bin/bash

I=1
while [ $I -le 10 ] 
	do
	for S in 0.01 0.1 0.5
	do
    	for Z in 0.01 0.1 0.5
    	do
    		for U in 0.1 0.5 0.9
    		do
      			for O in 0.1 0.5 0.9
      			do
					sbatch /scratch/axiom/FAC/FBM/DBC/nsalamin/default/tgaboria/ABM/Scripts/sim_ABM_best-S$S\-Z$Z\-U$U\-O$O\-I$I\.sh
      			done
    		done
    	done
  	done
  I=$(( $I + 1 ))
done