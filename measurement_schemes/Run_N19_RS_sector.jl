#!/bin/sh
# ########## Begin Slurm header ##########
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=50:00:00
#SBATCH --mem=350gb
#SBATCH --cpus-per-task=48
#SBATCH --job-name=otoc_simulation_N19_RS_sector
#SBATCH --output="otoc_simulation_N19_RS_sector-%j.out"
########### End Slurm header ##########
#=
# load modules
# not needed - julia installed locally

# export JULIA_DEPOT_PATH=$SCRATCH
export ON_CLUSTER=1
exec julia --color=no --threads=96 --startup-file=no "$LCDIR/LightCones/measurement_schemes/shared_krylov_sector.jl" 19 1 90 RS z $1
=#
