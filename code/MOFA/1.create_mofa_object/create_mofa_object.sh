#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=250000
#SBATCH -p fast
#SBATCH -N 1
##SBATCH --mem=100000
#SBATCH --mail-user=guillepalou4@gmail.com
##SBATCH --error=ewas_914_2055_ind_%j_%x.err
#SBATCH --output=ewas_914_2055_ind_%j_MOFA_object_groups.out
#SBATCH --job-name=MOFAobject

singularity exec -B /projects/regicor:/projects/regicor /projects/regicor/Guille/TFM/MOFA/MOFA2_docker.img Rscript /projects/regicor/Guille/TFM/MOFA/1.create_mofa_object/create_mofa_object_all_individuals.R
