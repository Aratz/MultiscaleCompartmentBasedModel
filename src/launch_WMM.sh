#!/bin/bash
#SBATCH -p core -n 1 -t 1:00:00 -A snic2019-8-227
module load gcc/10.1.0
python WMM_rackham.py ${1} ${2} ${3}
