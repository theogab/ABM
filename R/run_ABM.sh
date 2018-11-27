#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e ABM-error.txt
#BSUB –u theo.gaboriau@unil.ch
#BSUB –N
#BSUB –q dbc-long

module add R/latest
Rscript /scratch/cluster/monthly/tgaboria/ABM/ABM_tests.R $1 10000
module rm R/latest