module Utils

using SparseArrays, LinearAlgebra, Plots
using SpinSymmetry, JLD2
using KrylovKit
import Dates
using ..LightCones

export logmsg, path_prefix, save, save_with_pos, save_TR
export SimulationParams, SimulationParamsED
export logrange

struct SimulationParams
    N::Int
    SHOTS::Int
    RANDOM_STATES::Bool
    N_RANDOM_STATES::Int
    OBSERVABLE::String
    DISORDER_PARAM::Float64
end

struct SimulationParamsED
    N::Int
    SHOTS::Int
    OBSERVABLE::String
    DISORDER_PARAM::Float64
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

function save_with_pos(data,params,positiondata,jobid,datapath)
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
        file["positiondata"] = positiondata
    end
end

function save_TR(fids, otocs, params, jobid, datapath)
    dname = dirname(datapath)
    if !isdir(dname)
        logmsg("Save directory: $dname does not exist. Creating!")
        mkpath(dname)
    end
    logmsg("Saving file: $datapath")
    jldopen(datapath, "w") do file
        file["fidelities"] = fids
        file["otocs"] = otocs
        file["params"] = params
        file["jobid"] = jobid
    end    
end

function logrange(min_exp,max_exp,max)
    res = vcat([0],[i*10.0^j for j in min_exp:max_exp for i in 1:9])
    filter!(e->e<=max,res)
    return res
end

end #module