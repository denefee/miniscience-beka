#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=RT
#SBATCH --cpus-per-task=4
#SBATCH --job-name=mdcode
#SBATCH --comment="MD code run"

./mdlj -N 512 -ns 10000 -dt 0.001 -rho 0.8 -outf 100 -T0 1.0 -novelo -onefile