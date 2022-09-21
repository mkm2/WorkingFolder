#!/bin/sh
# ########## Begin Slurm header ##########
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=14:00:00
#SBATCH --mem=200gb
#SBATCH --cpus-per-task=48
#SBATCH --job-name=otoc_simulation_N15_RS_sector_log
#SBATCH --output="otoc_simulation_N15_RS_sector_log-%j.out"
########### End Slurm header ##########
#=
# load modules
# not needed - julia installed locally

# export JULIA_DEPOT_PATH=$SCRATCH
export ON_CLUSTER=1
exec julia --color=no --threads=96 --startup-file=no "$LCDIR/LightCones/measurement_schemes/shared_krylov_sector_log.jl" 15 1 1 RS z $1
=#
