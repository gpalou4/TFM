#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=125000
#SBATCH -p long
#SBATCH -N 1
##SBATCH --mem=100000
##SBATCH --mail-user=guillepalou4@gmail.com
##SBATCH --error=914_2055_ind_%j_%x.err
#SBATCH --output=914_2055_ind_%j_MOFA_nogroups_prediction.out
#SBATCH --job-name=MOFA.factor

singularity exec -B /projects/regicor:/projects/regicor /projects/regicor/Guille/TFM/MOFA/prediction.img Rscript /projects/regicor/Guille/TFM/MOFA/3.downstream_analysis/ewas_nogroups_914_2055/prediction/prediction.R
