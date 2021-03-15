#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=250000
#SBATCH -p fast
#SBATCH -N 1
##SBATCH --mem=100000
#SBATCH --mail-user=guillepalou4@gmail.com
##SBATCH --error=regicor_ewas_391_ind_%j_%x.err
#SBATCH --output=regicor_ewas_391_ind_%j_MOFA_training_groups.out
#SBATCH --job-name=MOFAtraining

singularity exec -B /projects/regicor:/projects/regicor /projects/regicor/Guille/TFM/MOFA/MOFA2_docker.img Rscript /projects/regicor/Guille/TFM/MOFA/2.training_model/regicor/training_model.R