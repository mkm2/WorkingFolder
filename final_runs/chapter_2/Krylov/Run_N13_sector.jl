#!/bin/sh
# ########## Begin Slurm header ##########
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10:00:00
#SBATCH --mem=250gb
#SBATCH --cpus-per-task=48
#SBATCH --job-name=otoc_simulation_N13_Krylov_sector
#SBATCH --output="otoc_simulation_N13_Krylov_sector-%j.out"
########### End Slurm header ##########
#=
# load modules
# not needed - julia installed locally

# export JULIA_DEPOT_PATH=$SCRATCH
export ON_CLUSTER=1
exec julia --color=no --threads=96 --startup-file=no "$LCDIR/WorkingFolder/final_runs/chapter_2/Krylov/shared_krylov_sector.jl" 13 1 1000 RS z $1
=#
