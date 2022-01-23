## environment
import Pkg
Pkg.activate(".")

if haskey(ENV, "ON_CLUSTER")
    @eval using MKL
    println("Added MKL.jl")
end

Pkg.instantiate(; io=stdout)
Pkg.status(; io=stdout)

import Dates
using SpinSymmetry, Random, LinearAlgebra
using LightCones
using Plots


#Location of LOGS
LOGS = get(ENV, "LOGS", "")

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
println("#BLAS config:               $(BLAS.get_config())")

@show ARGS

N = parse(Int, ARGS[1])
#SHOTS = parse(Int, ARGS[2]) #$(date '+%Y-%m-%d')

## constants and ARGS

#BLOCK = div(N-1,2)
#BASIS = SymmetrizedBasis(zbasis(N, BLOCK), [], [])

LOCATION = joinpath(LOGS,"LightCones",Dates.format(Dates.today(), "yyyy-mm-dd"))

@show LOCATION
@show N

logmsg("*"^10 * "Running simulation" * "*"^10)

H = xxz(N,6)
ψ0 = normalize!(ones(2^N))

op1 = single_spin_op(σz,5,N)
op2 = single_spin_op(σz,1,N)

trange = 0:0.1:5

corr = zeros(Float64,length(trange))
@time for (ti,t) in enumerate(trange)
    corr[ti] = 2-2*otoc(H, op1, op2, t, ψ0)
end
print(corr)

logmsg("*"^10 * "Saving" * "*"^10)
save(corr, joinpath(LOCATION,'test.jld2'))