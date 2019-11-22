#!/bin/bash -l

module add Bioinformatics/Software/vital-it 
module add R/latest
Rscript /scratch/axiom/FAC/FBM/DBC/nsalamin/default/tgaboria/ABM/Scripts/install.R
module rm R/latest