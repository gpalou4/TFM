#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=125000
#SBATCH -p long
#SBATCH -N 1
##SBATCH --mem=100000
##SBATCH --mail-user=guillepalou4@gmail.com
##SBATCH --error=cpgs_shared_%j_%x.err
#SBATCH --output=cpgs_shared_%j.out
#SBATCH --job-name=MOFA.analysis

singularity exec -B /projects/regicor:/projects/regicor /projects/regicor/Guille/TFM/MOFA/MOFA2_docker.img Rscript /projects/regicor/Guille/TFM/cpgs_shared.R
