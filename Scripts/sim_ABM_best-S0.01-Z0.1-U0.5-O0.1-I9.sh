#!/bin/bash 
#SBATCH --account=nsalamin_default 
#SBATCH --workdir=/scratch/axiom/FAC/FBM/DBC/nsalamin/default/tgaboria/ABM
#SBATCH --partition=ax-normal 
#SBATCH --error=err_sim_ABM_best-S0.01-Z0.1-U0.5-O0.1-I9.txt 
#SBATCH --output=out_sim_ABM_best-S0.01-Z0.1-U0.5-O0.1-I9.txt 
#SBATCH --time=4:00:00  
module add Bioinformatics/Software/vital-it 
module load R/latest 
Rscript Scripts/sim_tree_ABM_best.R 0.01 0.1 0.5 0.1 9