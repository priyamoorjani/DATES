#!/bin/sh
#SBATCH -N 1
#SBATCH --partition=short
#SBATCH --ntasks-per-node=1
#SBATCH -t 2:00:00 --mem 32000
#SBATCH -J post 
#SBATCH --mail-user=moorjani@gmail.com

/home/pm82/bin/convertf -p par.convertf

echo "End successfully"

