#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=125000
#SBATCH -p long
#SBATCH -N 1
##SBATCH --mem=100000
#SBATCH --mail-user=guillepalou4@gmail.com
##SBATCH --error=1110_ind_%j_%x.err
#SBATCH --output=1110_ind_%j_batch_effect.out
#SBATCH --job-name=BATCH

#singularity pull shub://regicor/mrna
#singularity pull --name mrna.sif shub://regicor/mrna0

singularity exec -B /projects/regicor:/projects/regicor /projects/regicor/Guille/TFM/transcriptomics/QC/gen_exp_image.sif Rscript /projects/regicor/Guille/TFM/transcriptomics/QC/6.batch_effect/batch_effect.R
