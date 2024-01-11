#!/bin/bash
#PBS -l select=1:ncpus=1:mpiprocs=1:mem=40G
#PBS -l walltime=24:00:00
#PBS -P 13003083
#PBS -q normal
#PBS -N count_CNV_per_window
#PBS -o log/
#PBS -j oe

# Enter where the job script is submitted.
cd $PBS_O_WORKDIR || exit $? 
# if no log folder, then make one
[ -d log ] || mkdir log  

# load R environment into the cluster
module load r/4.2.0

# Run R script.
R CMD BATCH count_cnvs_per_sliding_window.R
