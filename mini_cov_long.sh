#!/bin/bash

#SBATCH --job-name=artic
#SBATCH --account="" 
#SBATCH --cpus-per-task=48 
#SBATCH --nodes=1
#SBATCH --time=12:00:00 
#SBATCH --partition="ct560" # or ctest if time<0.5 hr
#SBATCH --mem=96G
#SBATCH -o %j.log           #path to std output
#SBATCH -e %j.err           #path to std err 

sh /home/ylllab2021/script/mini_cov.sh
