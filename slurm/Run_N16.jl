#!/bin/sh
# ########## Begin Slurm header ##########
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --mem=192gb
#SBATCH --cpus-per-task=48
#SBATCH --job-name=krylovkit-N16-a6
#SBATCH --output="$LOGS/LightCones/test/krylovkit-N_16-alpha_6-%j.out"
########### End Slurm header ##########
#=
# load modules
# not needed - julia installed locally

# export JULIA_DEPOT_PATH=$SCRATCH
export ON_CLUSTER=1
exec julia --color=no --threads=96 --startup-file=no "$LCDIR/LightCones/slurm/shared_krylov.jl" 16
=#