#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=250000
#SBATCH -p fast
#SBATCH -N 1
##SBATCH --mem=100000
#SBATCH --mail-user=guillepalou4@gmail.com
#SBATCH --error=50_%j_fast.err
#SBATCH --output=50_%j_fast.out
#SBATCH --job-name=input
##SBATCH -w, --nodelist=rslurmd-2.rslurmd

#singularity pull shub://regicor/mrna
#singularity pull --name mrna.sif shub://regicor/mrna

singularity exec -B /projects/regicor:/projects/regicor /projects/regicor/Guille/TFM/transcriptomics/QC/gen_exp_image.sif Rscript /projects/regicor/Guille/TFM/transcriptomics/QC/1.read_input/read_input.R
