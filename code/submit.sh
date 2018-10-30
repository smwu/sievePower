#!/bin/bash
#
#SBATCH --array=0-9
#SBATCH --job-name=markPower
#SBATCH --output=slurm_%a.out
export R_LIBS_USER=/home/tholzman/R/Library
ml R/3.3.3-foss-2016b-fh2
/app/easybuild/software/R/3.3.3-foss-2016b-fh2/lib/R/bin/Rscript --vanilla slurm_run.R
