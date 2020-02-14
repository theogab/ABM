#!/bin/bash 
#SBATCH --account=nsalamin_default 
#SBATCH --workdir=/scratch/axiom/FAC/FBM/DBC/nsalamin/default/tgaboria/ABM
#SBATCH --partition=ax-normal 
#SBATCH --error=err_abm_gibbs_sm-0.01_sx-0.1_nu-0.5_om-0.1_sim-5.txt 
#SBATCH --output=out_abm_gibbs_sm-0.01_sx-0.1_nu-0.5_om-0.1_sim-5.txt 
#SBATCH --time=4:00:00  
module add Bioinformatics/Software/vital-it 
module load R/latest 
Rscript Scripts/test_abm_gibbs.R 0.01 0.1 0.5 0.1 5