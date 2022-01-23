module Utils

using SparseArrays, LinearAlgebra, Plots
using SpinSymmetry, JLD2
using KrylovKit
import Dates
using ..LightCones

export logmsg, path_prefix, save

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

function save(data, datapath)
    dname = dirname(datapath)
    if !isdir(dname)
        logmsg("Save directory: $dname does not exist. Creating!")
        mkpath(dname)
    end
    logmsg("Saving file: $datapath")
    JLD2.jldsave(datapath; data) 
end

end #module