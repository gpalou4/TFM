#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=250000
#SBATCH -p fast
#SBATCH -N 1
##SBATCH --mem=100000
#SBATCH --mail-user=guillepalou4@gmail.com
##SBATCH --error=1200_ind_%j_%x.err
#SBATCH --output=1200_ind_%j_NUSEplots.out
#SBATCH --job-name=NUSEplots

#singularity pull shub://regicor/mrna
#singularity pull --name mrna.sif shub://regicor/mrna0

singularity exec -B /projects/regicor:/projects/regicor /projects/regicor/Guille/TFM/transcriptomics/QC/gen_exp_image.sif Rscript /projects/regicor/Guille/TFM/transcriptomics/QC/2.plots/NUSE_plots/plot_NUSE.R
