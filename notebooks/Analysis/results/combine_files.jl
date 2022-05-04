using LinearAlgebra
using JLD2
using LightCones

function combine_files(files,path,new_file)
	jobids = Vector{String}(undef,length(files))
	params = Vector{Any}(undef,length(files))
	data = Vector{Array{Float64,3}}(undef,length(files))
	for (i,f) in enumerate(files)
		jobids[i] = load(path*f,"jobid")
		params[i] = load(path*f,"params")
		data[i] = load(path*f,"data")
	end
	jldopen(path*new_file, "w") do file
        file["data"] = data
        file["params"] = params
        file["jobid"] = jobids
    end
end


path = pwd()*"/"
files = Vector{String}([x for x in ARGS[1:length(ARGS)-1]])
new_file = ARGS[length(ARGS)]

combine_files(files,path,new_file)
