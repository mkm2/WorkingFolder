#!/bin/sh
# ########## Begin Slurm header ##########
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=30:00:00
#SBATCH --mem=40gb
#SBATCH --cpus-per-task=48
#SBATCH --job-name=otoc_simulation_N14_short
#SBATCH --output="otoc_simulation_N14_short-%j.out"
########### End Slurm header ##########
#=
# load modules
# not needed - julia installed locally

# export JULIA_DEPOT_PATH=$SCRATCH
export ON_CLUSTER=1
exec julia --color=no --threads=96 --startup-file=no "$LCDIR/WorkingFolder/slurm/shared_krylov.jl" 14 1 0 z $1
=#
