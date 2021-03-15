#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=125000
#SBATCH -p long
#SBATCH -N 1
##SBATCH --mem=100000
##SBATCH --mail-user=guillepalou4@gmail.com
##SBATCH --error=ewas_914_2055_ind_%j_%x.err
#SBATCH --output=ewas_914_2055_ind_%j_MOFA_nogroups_cov_statistical_tests.out
#SBATCH --job-name=MOFA.analysis

singularity exec -B /projects/regicor:/projects/regicor /projects/regicor/Guille/TFM/MOFA/MOFA2_docker.img Rscript /projects/regicor/Guille/TFM/MOFA/3.downstream_analysis/ewas_nogroups_914_2055/statistical_tests/statistical_tests.R
