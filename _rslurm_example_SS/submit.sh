#!/bin/bash
#
#SBATCH --array=0-141
#SBATCH --cpus-per-task=1
#SBATCH --job-name=example_SS
#SBATCH --output=slurm_%a.out
#SBATCH --time=23:59:00
#SBATCH --mem=8192
#SBATCH --partition=short-serial
#SBATCH --error=%a.err
/usr/local/lib/R/bin/Rscript --vanilla slurm_run.R
