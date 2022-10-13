## environment
import Pkg
Pkg.activate(".")

using LinearAlgebra

if haskey(ENV, "ON_CLUSTER")
    @eval using MKL
    println("Added MKL.jl")
else
    BLAS.set_num_threads(1)
end

Pkg.instantiate(; io=stdout)
Pkg.status(; io=stdout)

import Dates
using SparseArrays, ThreadedSparseArrays
using SpinSymmetry, Random
using LightCones

## Environment
LOGS = get(ENV, "LOGS", "")
JOBID = get(ENV, "SLURM_JOB_ID", "")

logmsg("*"^10*"RANDOM FIELDS"*"*"^10)
println("shared_krylov_sector.jl")

println("Working Directory:          $(pwd())" )
println("SLURM Directory:            $(get(ENV, "SLURM_SUBMIT_DIR", "")) ")
println("Running on host:            $(gethostname())" )
println("Job id:                     $(get(ENV, "SLURM_JOB_ID", ""))" )
println("Job name:                   $(get(ENV, "SLURM_JOB_NAME", ""))" )
println("Number of nodes allocated:  $(get(ENV, "SLURM_JOB_NUM_MODES", ""))" )
println("Number of cores allocated:  $(get(ENV, "SLURM_NTASKS", ""))" )
println("#threads of Julia:          $(Threads.nthreads())")
println("#threads of BLAS:           $(BLAS.get_num_threads())")
#println("#BLAS config:               $(BLAS.get_config())")

@show ARGS

## Constants and ARGS

N = parse(Int, ARGS[1])
SHOTS = parse(Int, ARGS[2])
OBSERVABLE = ARGS[3]
DISORDER_PARAM = parse(Float64, ARGS[4])

LOCATION = joinpath(LOGS,"LightCones",Dates.format(Dates.today(), "yyyy-mm-dd"))

@show LOCATION
@show N
@show SHOTS
@show OBSERVABLE
@show DISORDER_PARAM

params = SimulationParamsED(N,SHOTS,OBSERVABLE,DISORDER_PARAM)

logmsg("*"^10 * "Running simulation" * "*"^10)

#Set up simulation parameters

s = 10
trange = 10. .^LinRange(-3,8,100)
logmsg("trange = ",trange)

i = div(N,2)+1
k = div(N-1,2)+1 #largest sector
d = basissize(symmetrized_basis(N,k))

A = convert(SparseMatrixCSC{ComplexF64,Int64},symmetrize_operator(single_spin_op(σz,i,N),N,k))
logmsg("A=σz")

if OBSERVABLE == "x"
    B = σx
elseif OBSERVABLE == "y"
    B = σy
elseif OBSERVABLE == "z"
    B = σz
end
B = convert(SparseMatrixCSC{ComplexF64,Int64},B)

H = symmetrize_operator(xxz(N,6),N,k)
logmsg(typeof(H))

#Start simulation

otocs = zeros(length(trange),N,SHOTS)
H_tot = Vector{SparseMatrixCSC{ComplexF64,Int64}}([spzeros(d,d) for l in 1:SHOTS])
Threads.@threads for shot in 1:SHOTS
    H_tot[shot] = H + field_term(DISORDER_PARAM,N,k)
    logmsg("Created Hamiltonian for Shot $(shot)")
    otocs[:,:,shot] = Diag_OTOC(Matrix(H_tot[shot]),A,B,trange,N,s,k,d)
    logmsg("Completed Shot $(shot)")
end

logmsg("*"^10*"Simulation completed!"*"*"^10)

logmsg("*"^10 * "Saving" * "*"^10)
save(otocs, params, JOBID, joinpath(LOCATION,"$(JOBID)_N$(N)_ED.jld2"))
logmsg("*"^10 * "Run completed!" * "*"^10)
