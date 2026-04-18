#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-100
#SBATCH --mail-type=NONE
#SBATCH --mail-user=aron0064@umn.edu
#SBATCH -A mfiecas
#SBATCH -o /panfs/jay/groups/29/mfiecas/aron0064/bkmrLoD/LogFiles/%A_%a.out
#SBATCH -e /panfs/jay/groups/29/mfiecas/aron0064/bkmrLoD/LogFiles/%A_%a.err
n=$1
lod_quantile=$2
exposure_dist=$3
h_dist=$4
date
path="$HOME/bkmrLoDsim"
cd $path/Routputs
start_time=$(date +%s)
module load R/4.4.0
Rscript $path/Rcode/bmkrAug.R $n $lod_quantile $exposure_dist $h_dist
finish_time=$(date +%s)
elapsed_time=$((finish_time  - start_time))

((sec=elapsed_time%60, elapsed_time/=60, min=elapsed_time%60, hrs=elapsed_time/60))
timestamp=$(printf "Total time taken - %d hours, %d minutes, and %d seconds." $hrs $min $sec)
echo $timestamp

