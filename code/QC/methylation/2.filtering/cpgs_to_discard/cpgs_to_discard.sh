#!/bin/bash

#SBATCH -n 1 
#SBATCH --cpus-per-task=4
#SBATCH -p fast
#SBATCH -N 1
#SBATCH --mem=100000

singularity exec -B /projects/regicor:/projects/regicor /projects/regicor/METIAM/EPIC/QC/r_ewas6.img Rscript /projects/regicor/Guille/TFM/methylation/QC/2.filtering/cpgs_to_discard/cpgs_to_discard.R
