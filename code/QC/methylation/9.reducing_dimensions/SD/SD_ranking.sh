#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=125000
#SBATCH -p long
#SBATCH -N 1
##SBATCH --mem=100000

singularity exec -B /projects/regicor:/projects/regicor /projects/regicor/METIAM/EPIC/QC/r_ewas6.img Rscript /projects/regicor/Guille/TFM/methylation/QC/10.reducing_dimensions/SD/SD_ranking.R
