#!/bin/bash
#SBATCH -p node -n 1 -t 59:00 -A snic2019-8-227
module load gcc/9.2

python run_inference.py ${1} ${2}
