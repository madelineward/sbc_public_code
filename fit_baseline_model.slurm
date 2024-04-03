#!/bin/bash
####### Reserve computing resources #############
#SBATCH --job-name=base_models.run
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=8:30:00
#SBATCH --mem 4G
#SBATCH --partition=cpu2019,cpu2021-bf24,cpu2021,cpu2022,cpu2023,cpu2022-bf24
#SBATCH --array=1-20
#SBATCH --output=./out/base.%A_%a.out
#SBATCH --error=./err/base.%A_%a.error

####### Run your script #########################
/home/madeline.ward1/julia-1.6.3/bin/julia fit_baseline_model.jl