#!/bin/sh
# ########## Begin Slurm header ##########
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=05:00:00
#SBATCH --mem=192gb
#SBATCH --cpus-per-task=2
#SBATCH --job-name=otoc_simulation_N16
#SBATCH --output="otoc_simulation_N16-%j.out"
########### End Slurm header ##########
#=
# load modules
# not needed - julia installed locally

# export JULIA_DEPOT_PATH=$SCRATCH
export ON_CLUSTER=1
exec julia --color=no --threads=4 --startup-file=no "$LCDIR/LightCones/slurm/shared_krylov.jl" 16 10 0 z $1
=#
