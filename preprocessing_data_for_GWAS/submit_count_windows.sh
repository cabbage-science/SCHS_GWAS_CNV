#!/bin/bash
#PBS -l select=1:ncpus=1:mpiprocs=1:mem=40G
#PBS -l walltime=24:00:00
#PBS -P 13003083
#PBS -q normal
#PBS -N count_CNV_per_window

# load R environment into the cluster
module load r/4.2.0

# Run R script.
R CMD BATCH count_cnvs_per_sliding_window.R
