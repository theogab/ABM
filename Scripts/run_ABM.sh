#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e ABM-error.txt
#BSUB –u theo.gaboriau@unil.ch
#BSUB –N
#BSUB –q dbc-long

module add R/latest
Rscript /scratch/cluster/monthly/tgaboria/ABM/R/sim_tree_ABM.R $1 $2 0.002 100
Rscript /scratch/cluster/monthly/tgaboria/ABM/R/sim_tree_ABM.R $1 $2 0.02 100
Rscript /scratch/cluster/monthly/tgaboria/ABM/R/sim_tree_ABM.R $1 $2 0.2 100
module rm R/latest