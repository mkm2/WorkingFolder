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
using SimLib, XXZNumerics


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
#SHOTS = parse(Int, ARGS[2])

## constants and ARGS

#BLOCK = div(N-1,2)
#BASIS = SymmetrizedBasis(zbasis(N, BLOCK), [], [])

LOCATION = joinpath(LOGS,"LightCones",Dates.format(Dates.today(), "yyyy-mm-dd"))

@show LOCATION
@show N

H = xxz()

model = RandomPositionsXXZWithXField(pdd, PowerLaw(ALPHA), [0.0], nothing, BASIS)

edd = EDDataDescriptor(model, DIAGTYPE, LOCATION)

logmsg("*"^10 * "Saving" * "*"^10)
save.(edata)