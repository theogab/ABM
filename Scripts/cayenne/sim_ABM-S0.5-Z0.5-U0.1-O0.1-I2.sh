#!/bin/bash 
#SBATCH --account=nsalamin_default 
#SBATCH --workdir=/scratch/axiom/FAC/FBM/DBC/nsalamin/default/tgaboria/ABM
#SBATCH --partition=ax-normal 
#SBATCH --error=err_sim_ABM-S0.5-Z0.5-U0.1-O0.1-I2.txt 
#SBATCH --output=out_sim_ABM-S0.5-Z0.5-U0.1-O0.1-I2.txt 
#SBATCH --time=4:00:00  
module add Bioinformatics/Software/vital-it 
module load R/latest 
Rscript Scripts/sim_tree_ABM.R 0.5 0.5 0.1 0.1 2