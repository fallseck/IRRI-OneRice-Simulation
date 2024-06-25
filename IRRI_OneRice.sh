#!/bin/bash
#SBATCH -J IRRI-OneRice #Job name
#SBATCH -p agap_normal #Choose a Partition to use 
#SBATCH --cpus-per-task=2 #cpus per task 5
#SBATCH --mem-per-cpu=8GB #Memory size
#SBATCH -a 1-80
#SBATCH --mail-type=BEGIN,END,FAIL #Notifications 
#SBATCH --mail-user=f.seck@irri.org #Address to send
#SBATCH -t 47:59:59 #Time limitation

module load R/packages/4.1.0

Rscript master.R
