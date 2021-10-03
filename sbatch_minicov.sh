#!/bin/bash

#SBATCH --job-name=artic
#SBATCH --account="GOV110078" 
#SBATCH --cpus-per-task=48 
#SBATCH --nodes=1
#SBATCH --time=48:00:00 
#SBATCH --partition="ct56" # or ctest if time<0.5 hr
#SBATCH --mem=96G
#SBATCH -o %j.log           #path to std output
#SBATCH -e %j.err           #path to std err 

sh mini_cov.sh
