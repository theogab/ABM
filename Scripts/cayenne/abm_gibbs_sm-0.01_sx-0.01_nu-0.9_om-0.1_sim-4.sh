#!/bin/bash 
#SBATCH --account=nsalamin_default 
#SBATCH --workdir=/scratch/axiom/FAC/FBM/DBC/nsalamin/default/tgaboria/ABM
#SBATCH --partition=ax-normal 
#SBATCH --error=err_abm_gibbs_sm-0.01_sx-0.01_nu-0.9_om-0.1_sim-4.txt 
#SBATCH --output=out_abm_gibbs_sm-0.01_sx-0.01_nu-0.9_om-0.1_sim-4.txt 
#SBATCH --time=4:00:00  
module add Bioinformatics/Software/vital-it 
module load R/latest 
Rscript Scripts/test_abm_gibbs.R 0.01 0.01 0.9 0.1 4