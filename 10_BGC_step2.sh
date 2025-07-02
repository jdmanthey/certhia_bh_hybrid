#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=bgc
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=4
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-100

source activate bgc

# set the number of SNPs to run per job
# can figure out in R with: ceiling(n_snps / n_array_jobs)
PER_TASK=371

# calculate the starting and ending SNPs for this task based
# on the SLURM task and the number of SNPs per task
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * PER_TASK ))

# for the last job, change the END_NUM variable to the last SNP
if (($SLURM_ARRAY_TASK_ID == 100 )); then END_NUM=37017; fi

# print the task and range of SNPs for this run
echo This is SLURM task $SLURM_ARRAY_TASK_ID, which will do BGC estimation for SNPs $START_NUM to $END_NUM

# run the R script to fit clines for each SNP
Rscript _BGC_step2_helper.r $START_NUM $END_NUM $SLURM_ARRAY_TASK_ID




