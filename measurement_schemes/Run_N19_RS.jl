#!/bin/sh
# ########## Begin Slurm header ##########
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=100gb
#SBATCH --cpus-per-task=48
#SBATCH --job-name=otoc_simulation_N19_RS
#SBATCH --output="otoc_simulation_N19_RS-%j.out"
########### End Slurm header ##########
#=
# load modules
# not needed - julia installed locally

# export JULIA_DEPOT_PATH=$SCRATCH
export ON_CLUSTER=1
exec julia --color=no --threads=96 --startup-file=no "$LCDIR/WorkingFolder/measurement_schemes/shared_krylov.jl" 19 1 10 RS z $1
=#
