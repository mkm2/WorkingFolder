#!/bin/sh
# ########## Begin Slurm header ##########
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=250gb
#SBATCH --cpus-per-task=48
#SBATCH --job-name=otoc_simulation_N12_ED_tfim
#SBATCH --output="otoc_simulation_N12_ED_tfim-%j.out"
########### End Slurm header ##########
#=
# load modules
# not needed - julia installed locally

# export JULIA_DEPOT_PATH=$SCRATCH
export ON_CLUSTER=1
exec julia --color=no --threads=96 --startup-file=no "$LCDIR/WorkingFolder/Strasbourg/shared_ED_tfim.jl" 12 100 z $1
=#
