#!/bin/bash
#SBATCH -p core -n 1 -t 12:00:00 -A snic2019-8-227
python smoldyn_rackham.py ${1} ${2} ${3}
