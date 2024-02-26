#!/bin/bash -l
#SBATCH -c 64
#SBATCH --time=0-00:30:00
#SBATCH -p batch

module load lang/Julia/1.8
module load math/Gurobi

julia example_calibrate_para_n_imp.jl
