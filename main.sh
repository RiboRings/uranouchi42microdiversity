#!/bin/bash

#PBS -q cdb
#PBS -o stdout.txt
#PBS -e err.txt
#PBS -N Rscript
#PBS -l mem=200gb
#PBS -l ncpus=15
#PBS -l nice=5
#PBS -m ea
#PBS -M giuliobene2000@gmail.com

source /etc/profile.d/modules.sh
cd /lustre1/aptmp/giulio/uranouchi42microdiversity

echo "loading module"
module load R
echo "module loaded"

Rscript main.R

