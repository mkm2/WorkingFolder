#!/bin/sh
# ########## Begin Slurm header ##########
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=7-0
#SBATCH --mem=350gb
#SBATCH --cpus-per-task=48
#SBATCH --job-name=otoc_simulation_N17_EDTR_sector
#SBATCH --output="otoc_simulation_N17_EDTR_sector-%j.out"
########### End Slurm header ##########
#=
# load modules
# not needed - julia installed locally

# export JULIA_DEPOT_PATH=$SCRATCH
export ON_CLUSTER=1
exec julia --color=no --threads=96 --startup-file=no "$LCDIR/WorkingFolder/final_runs/chapter_2/ED_tr/shared_krylov_sector.jl" 17 1 z $1
=#
