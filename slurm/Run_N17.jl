#!/bin/sh
# ########## Begin Slurm header ##########
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=100:00:00
#SBATCH --mem=100gb
#SBATCH --cpus-per-task=48
#SBATCH --job-name=otoc_simulation_N17_short
#SBATCH --output="otoc_simulation_N17_short-%j.out"
########### End Slurm header ##########
#=
# load modules
# not needed - julia installed locally

# export JULIA_DEPOT_PATH=$SCRATCH
export ON_CLUSTER=1
exec julia --color=no --threads=96 --startup-file=no "$LCDIR/LightCones/slurm/shared_krylov.jl" 17 1 0 z $1
=#
