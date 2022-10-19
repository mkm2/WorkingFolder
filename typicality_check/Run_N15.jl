#!/bin/sh
# ########## Begin Slurm header ##########
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=180gb
#SBATCH --cpus-per-task=48
#SBATCH --job-name=otoc_simulation_N15_RPS
#SBATCH --output="otoc_simulation_N15_RPS-%j.out"
########### End Slurm header ##########
#=
# load modules
# not needed - julia installed locally

# export JULIA_DEPOT_PATH=$SCRATCH
export ON_CLUSTER=1
exec julia --color=no --threads=96 --startup-file=no "$LCDIR/WorkingFolder/typicality_check/shared_krylov.jl" 15 1 100 z $1
=#
