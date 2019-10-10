#!/bin/bash
#SBATCH --partition=cpu_short
#SBATCH --job-name=abcHU
#SBATCH --mem-per-cpu=9000
#SBATCH --time=11:50:00
#SBATCH --tasks=4
#SBATCH --cpus-per-task=1
#SBATCH --nodes=2
#SBATCH --output=Rdata


nohup Rscript $*

