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
SEQ_REPS = parse(Int, ARGS[7])
SEQ_NAME = ARGS[8]

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

tmax = 0.5
#trange = 10. .^LinRange(-3,1,100)
trange = 0:0.1:20
logmsg("trange = ",trange)

i = div(N,2)+1
A = convert(SparseMatrixCSC{ComplexF64,Int64},single_spin_op(σz,i,N))
logmsg("A = σz")

if OBSERVABLE == "x"
    B = σx
elseif OBSERVABLE == "y"
    B = σy
elseif OBSERVABLE == "z"
    B = σz
end
B = convert(SparseMatrixCSC{ComplexF64,Int64},B)

H = xxz_pbc(N,6)
logmsg("alpha=6")
logmsg("PBC")

if MULT_RANDOM_STATES == false
    if TYPE_OF_RS == "RS"
        ψ0 = random_state(N)#normalize!(ones(2^N))
    elseif TYPE_OF_RS == "RPS"
        ψ0 = random_product_state(N)
    elseif TYPE_OF_RS == "BS"
        ψ0 = random_bitstring_state(N)
    else
        throw("No states sampled. Wrong input.")
    end
    logmsg("Sampled 1 random initial state")
else
    ψs = zeros(ComplexF64,2^N,N_RANDOM_STATES)
    if TYPE_OF_RS == "RS"
        for s in 1:N_RANDOM_STATES
            ψs[:,s] = random_state(N)
        end
    elseif TYPE_OF_RS == "RPS"
        for s in 1:N_RANDOM_STATES
            ψs[:,s] = random_product_state(N)
        end
    elseif TYPE_OF_RS == "BS"
        for s in 1:N_RANDOM_STATES
            ψs[:,s] = random_bitstring_state(N,OBSERVABLE)
        end
    else
        throw("No states sampled. Wrong input.")
    end
end

#Start simulation
if MULT_RANDOM_STATES == false
    signs = signs_of_eigenstate(B,ψ0,N)
    fidelities = zeros(length(trange),SHOTS)
    oto_commutators = zeros(length(trange),N,SHOTS)
    H_tot = Vector{SparseMatrixCSC{Float64,Int64}}([spzeros(2^N,2^N) for l in 1:SHOTS])
    #H_tot = Vector{Adjoint{Float64, ThreadedSparseMatrixCSC{Float64, Int64, SparseMatrixCSC{Float64, Int64}}}}([ThreadedSparseMatrixCSC(spzeros(2^N,2^N))' for l in 1:4])
    Threads.@threads for shot in 1:SHOTS
        H_tot[shot] = H + field_term(DISORDER_PARAM,N)
        logmsg("Created Hamiltonian for Shot $(shot)")
        #H_tot[shot] = ThreadedSparseMatrixCSC(H + field_term(DISORDER_PARAM,N))'
        @time fidelities[:,shot] = fidelity_timereversal(H_tot[shot],trange,ψ0,SEQ_NAME,SEQ_REPS,N,"Krylov",tmax)
        @time oto_commutators[:,:,shot] = oto_commutator_timereversal(H_tot[shot],A,B,trange,ψ0,SEQ_NAME,SEQ_REPS,N,"Krylov",tmax)
        logmsg("Completed Shot $(shot)")
    end
else
    fidelities = zeros(length(trange),SHOTS,N_RANDOM_STATES)
    oto_commutators = zeros(length(trange),N,SHOTS,N_RANDOM_STATES)
    print("Using multiple states.\n")
    for shot in 1:SHOTS
        #H_tot[shot] = ThreadedSparseMatrixCSC(H + field_term(DISORDER_PARAM,N))'
        H_tot = H + field_term(DISORDER_PARAM,N)
        Threads.@threads for s in 1:N_RANDOM_STATES
            @time fidelities[:,shot,s] = fidelity_timereversal(H_tot,trange,ψs[:,s],SEQ_NAME,SEQ_REPS,N,"Krylov",tmax)
            @time oto_commutators[:,:,shot,s] = oto_commutator_timereversal(H_tot,A,B,trange,ψs[:,s],SEQ_NAME,SEQ_REPS,N,"Krylov",tmax)
            logmsg("Completed Shot $(shot), state $(s)")
        end
    end
end

logmsg("*"^10*"Simulation completed!"*"*"^10)

logmsg("*"^10 * "Saving" * "*"^10)
save_TR(fidelities, oto_commutators, params, JOBID, joinpath(LOCATION,"$(JOBID)_N$(N)_$(TYPE_OF_RS).jld2"))
logmsg("*"^10 * "Run completed!" * "*"^10)
