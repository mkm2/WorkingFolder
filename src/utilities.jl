module Utils

using SparseArrays, LinearAlgebra, Plots
using SpinSymmetry, JLD2
using KrylovKit
import Dates
using ..LightCones

export logmsg, path_prefix, save
export Maybe, SimulationParams

const Maybe{T} = Union{Missing, T} where T

struct SimulationParams
    N::Int
    SHOTS::Int
    RANDOM_STATES::Bool
    OBSERVABLE::String
    DISORDER_PARAM::Float64
    N_RANDOM_STATES::Maybe{Int}
end

function logmsg(msg...; doflush=false)
    println("[",Dates.now(), "]", msg...)
    doflush && flush(stdout)
end

function path_prefix(workspace)
    try
        abspath(joinpath(readchomp(`ws_find $workspace`)))
    catch e
        path_prefix()
    end
end

function save(data, params, jobid, datapath)
    dname = dirname(datapath)
    if !isdir(dname)
        logmsg("Save directory: $dname does not exist. Creating!")
        mkpath(dname)
    end
    logmsg("Saving file: $datapath")
    jldopen(datapath, "w") do file
        file["data"] = data
        file["params"] = params
        file["jobid"] = jobid
    end
end

end #module