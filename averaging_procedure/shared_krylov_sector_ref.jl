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
println("shared_krylov.jl")

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
N_RANDOM_STATES = parse(Int, ARGS[3])
if N_RANDOM_STATES == 0
    MULT_RANDOM_STATES = false
else
    MULT_RANDOM_STATES = true
end
TYPE_OF_RS = ARGS[4]
OBSERVABLE = ARGS[5]
DISORDER_PARAM = parse(Float64, ARGS[6])

#SHOTS = parse(Int, ARGS[2]) #$(date '+%Y-%m-%d')

#BLOCK = div(N-1,2)
#BASIS = SymmetrizedBasis(zbasis(N, BLOCK), [], [])

LOCATION = joinpath(LOGS,"LightCones",Dates.format(Dates.today(), "yyyy-mm-dd"))

@show LOCATION
@show N
@show SHOTS
@show MULT_RANDOM_STATES
@show N_RANDOM_STATES
@show TYPE_OF_RS
@show OBSERVABLE
@show DISORDER_PARAM

params = SimulationParams(N,SHOTS,MULT_RANDOM_STATES,N_RANDOM_STATES,OBSERVABLE,DISORDER_PARAM)

logmsg("*"^10 * "Running simulation" * "*"^10)

#Set up simulation parameters

δt = 0.1
tmax = 0.5
T = 5
#trange = logrange(-5,0,2)
trange = 0:δt:T
logmsg("trange = ",trange)

i = div(N,2)+1
k = div(N-1,2)+1 #largest sector
d = basissize(symmetrized_basis(N,k))

A = symmetrize_operator(single_spin_op(σz,i,N),N,k)
logmsg("A=σz")

if OBSERVABLE == "x"
    B = σx
elseif OBSERVABLE == "y"
    B = σy
elseif OBSERVABLE == "z"
    B = σz
end

print("before H\n")
H = symmetrize_operator(xxz(N,6),N,k)
print("after H\n")
print(string("memory allocated: ",Base.summarysize(H)))
print(string("memory of full H: ",Base.summarysize(xxz(N,6))))

#Start simulation
otocs = zeros(length(trange),N,SHOTS,N_RANDOM_STATES)
H_tot = Vector{SparseMatrixCSC{Float64,Int64}}([spzeros(d,d) for l in 1:SHOTS])
print("test\n")
for shot in 1:SHOTS
    #H_tot[shot] = ThreadedSparseMatrixCSC(H + field_term(DISORDER_PARAM,N))'
    H_tot[shot] = H + field_term(DISORDER_PARAM,N,k)
    Threads.@threads for s in 1:N_RANDOM_STATES
        @time otocs[:,:,shot,s] = otoc_spat(H_tot[shot],A,B,trange,random_state(N,d),N,k,tmax)
        logmsg("Completed Shot $(shot), state $(s)")
    end
end

logmsg("*"^10*"Simulation completed!"*"*"^10)

logmsg("*"^10 * "Saving" * "*"^10)
save(otocs, params, JOBID, joinpath(LOCATION,"$(JOBID)_N$(N)_$(TYPE_OF_RS).jld2"))
logmsg("*"^10 * "Run completed!" * "*"^10)
