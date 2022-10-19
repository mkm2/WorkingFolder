#!/bin/sh
# ########## Begin Slurm header ##########
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=250gb
#SBATCH --cpus-per-task=48
#SBATCH --job-name=otoc_simulation_N13_ED_nn
#SBATCH --output="otoc_simulation_N13_ED_nn-%j.out"
########### End Slurm header ##########
#=
# load modules
# not needed - julia installed locally

# export JULIA_DEPOT_PATH=$SCRATCH
export ON_CLUSTER=1
exec julia --color=no --threads=96 --startup-file=no "$LCDIR/WorkingFolder/exactdiag/shared_krylov_nn.jl" 13 50 z $1
=#
