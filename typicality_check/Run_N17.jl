#!/bin/sh
# ########## Begin Slurm header ##########
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=80gb
#SBATCH --cpus-per-task=10
#SBATCH --job-name=otoc_simulation_N17
#SBATCH --output="otoc_simulation_N17-%j.out"
########### End Slurm header ##########
#=
# load modules
# not needed - julia installed locally

# export JULIA_DEPOT_PATH=$SCRATCH
export ON_CLUSTER=1
exec julia --color=no --threads=20 --startup-file=no "$LCDIR/LightCones/typicality_check/shared_krylov.jl" 17 1 10 z $1
=#
