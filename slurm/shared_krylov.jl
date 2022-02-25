## environment
import Pkg
Pkg.activate(".")

if haskey(ENV, "ON_CLUSTER")
    @eval using MKL
    println("Added MKL.jl")
else
    BLAS.set_num_threads(1)
end

Pkg.instantiate(; io=stdout)
Pkg.status(; io=stdout)

import Dates
using SpinSymmetry, Random, LinearAlgebra
using LightCones

## Environment
LOGS = get(ENV, "LOGS", "")
JOBID = get(ENV, "SLURM_JOB_ID", "")

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
    RANDOM_STATES = false
else
    RANDOM_STATES = true
end
OBSERVABLE = ARGS[4]
DISORDER_PARAM = parse(Int, ARGS[5])

#SHOTS = parse(Int, ARGS[2]) #$(date '+%Y-%m-%d')

#BLOCK = div(N-1,2)
#BASIS = SymmetrizedBasis(zbasis(N, BLOCK), [], [])

LOCATION = joinpath(LOGS,"LightCones",Dates.format(Dates.today(), "yyyy-mm-dd"))

@show LOCATION
@show N
@show SHOTS
@show RANDOM_STATES
@show N_RANDOM_STATES
@show OBSERVABLE

if RANDOM_STATES == false
    params = SimulationParams(N,SHOTS,RANDOM_STATES,OBSERVABLE,DISORDER_PARAM)
else
    params = SimulationParams(N,SHOTS,RANDOM_STATES,N_RANDOM_STATES,OBSERVABLE,DISORDER_PARAM) 
end

logmsg("*"^10 * "Running simulation" * "*"^10)

#Set up simulation parameters

δt = 0.1
T = 5
trange = 0:δt:T

i = 3
A = single_spin_op(σz,i,N)

if OBSERVABLE == 'x'
    B = σx
else if OBSERVABLE == 'y'
    B = σy
else if OBSERVABLE == 'z'
    B = σz
end

H = xxz(N,6)
if RANDOM_STATES == false
    ψ0 = normalize!(ones(2^N))
else
    ψs = zeros(N_RANDOM_STATES,2^N)
    for s in 1:N_RANDOM_STATES
        ψs[s] = random_state(N)
    end
end

#Start simulation

if RANDOM_STATES == false
    otocs = zeros(length(trange),N,SHOTS)
    H_tot = spzeros(2^N,2^N,SHOTS)
    Threads.@threads for shot in 1:SHOTS
        H_tot[:,:,shot] = H + field_term(DISORDER_PARAM,N)
        otocs[:,:,shot] = otoc_spat(H_tot[shot],A,B,trange,ψ0,N,δt)
    end
else
    otocs = zeros(SHOTS,N_RANDOM_STATES,length(trange),N)
    H_tot = spzeros(2^N,2^N,SHOTS)
    Threads.@threads for shot in 1:SHOTS
        Threads.@threads for s in 1:N_RANDOM_STATES
            H_tot[:,:,shot] = H + field_term(DISORDER_PARAM,N)
            otocs[:,:,shot,s] = otoc_spat(H_tot[shot],A,B,trange,ψs[s],N,δt)
        end
    end
end

logmsg("*"^10*"Simulation completed!"*"*"^10)

logmsg("*"^10 * "Saving" * "*"^10)
save(otocs, params, JOBID, joinpath(LOCATION,"$(JOBID)_N$(N).jld2"))
logmsg("*"^10 * "Run completed!" * "*"^10)
