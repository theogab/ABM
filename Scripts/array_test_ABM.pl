#!/usr/bin/env perl

use File::Find;
$I=1;
while ($I<=10){
  foreach my $S (0.01, 0.1, 0.5){
    foreach my $Z (0.01, 0.1, 0.5){
    	foreach my $U (0.1, 0.5, 0.9){
      		foreach my $O (0.1, 0.5, 0.9){
         		my $handle = undef;
  		 		my $filename = "Scripts/sim_ABM_best-S$S\-Z$Z\-U$U\-O$O\-I$I\.sh";
  		 		my $encoding = ":encoding(UTF-8)";
  		 		open($handle, "> $encoding", $filename);
  		 		print {$handle} "#!/bin/bash \n#SBATCH --account=nsalamin_default \n#SBATCH --workdir=/scratch/axiom/FAC/FBM/DBC/nsalamin/default/tgaboria/ABM\n#SBATCH --partition=ax-normal \n#SBATCH --error=err_sim_ABM_best-S$S\-Z$Z\-U$U\-O$O\-I$I\.txt \n#SBATCH --output=out_sim_ABM_best-S$S\-Z$Z\-U$U\-O$O\-I$I\.txt \n#SBATCH --time=4:00:00  \nmodule add Bioinformatics/Software/vital-it \nmodule load R/latest \nRscript Scripts/sim_tree_ABM_best.R $S $Z $U $O $I";
  		 		close($handle);
      		}
    	}
    }
  }
  $I++;
}
exit(0);