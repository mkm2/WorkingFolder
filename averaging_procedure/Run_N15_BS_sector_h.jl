#!/bin/sh
# ########## Begin Slurm header ##########
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=30:00:00
#SBATCH --mem=180gb
#SBATCH --cpus-per-task=48
#SBATCH --job-name=otoc_simulation_N13_BS_sector_AVP_h
#SBATCH --output="otoc_simulation_N13_BS_sector_AVP_h-%j.out"
########### End Slurm header ##########
#=
# load modules
# not needed - julia installed locally

# export JULIA_DEPOT_PATH=$SCRATCH
export ON_CLUSTER=1
exec julia --color=no --threads=96 --startup-file=no "$LCDIR/WorkingFolder/averaging_procedure/shared_krylov_sector_h.jl" 17 50 $1 BS z $2
=#
