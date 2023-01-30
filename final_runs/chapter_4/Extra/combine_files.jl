### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ af1f73b4-268c-43d9-84f7-9f8d7090cdce
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 849b46da-e52e-4a3a-ba44-f5453cc06529
using LinearAlgebra,Plots,JLD2,Statistics,PlutoUI

# ╔═╡ a51014ea-2397-40a2-844e-11a5951e8a51
TableOfContents()

# ╔═╡ baf8383d-9370-4723-8418-c736cf0e04e4
function combine_files(files,path,new_file)
	jobids = Vector{String}(undef,length(files))
	params = Vector{Any}(undef,length(files))
	data = Vector{Array{Float64,3}}(undef,length(files))
	#positiondata = Vector{PositionData}(undef, length(files))
	for (i,f) in enumerate(files)
		jobids[i] = load(path*f,"jobid")
		params[i] = load(path*f,"params")
		#positiondata[i] = load(path*f,"positiondata")
		data[i] = load(path*f,"data")[:,:,:]
	end
	jldopen(path*new_file, "w") do file
        file["data"] = data
		#file["positiondata"] = positiondata
        file["params"] = params
        file["jobid"] = jobids
    end
end

# ╔═╡ 5b1b59eb-d749-45d2-9d04-19ea22fc3e22
names = ["nn","pl","k7","k10"]

# ╔═╡ 33a94d93-15a6-4a93-82e3-f1d21ad2df6b
path = pwd()

# ╔═╡ 2a8b7e95-ed80-422d-b655-9aad01353ec6
folders = ["/AL/nn/","/AL/pl/","/sector/k = 7/","/sector/k = 10/"]

# ╔═╡ 1e4d8b70-c322-4d79-82b7-a6209619ff88
begin
	for (i,folder) in enumerate(folders)
		name = names[i]
		print("---"*folder*"---\n\n")
			l = []
			for s in readdir(path*folders[i])
				if !occursin("otoc",s) && !occursin("combined",s)
					append!(l,[s])
				end
			end
			print("name $name:\n")
			print(l)
			print("\n\n")
			combine_files(l,path*folder,"combined_file_$name.jld2")
		end
end

# ╔═╡ Cell order:
# ╠═849b46da-e52e-4a3a-ba44-f5453cc06529
# ╠═af1f73b4-268c-43d9-84f7-9f8d7090cdce
# ╠═a51014ea-2397-40a2-844e-11a5951e8a51
# ╠═baf8383d-9370-4723-8418-c736cf0e04e4
# ╠═5b1b59eb-d749-45d2-9d04-19ea22fc3e22
# ╠═33a94d93-15a6-4a93-82e3-f1d21ad2df6b
# ╠═2a8b7e95-ed80-422d-b655-9aad01353ec6
# ╠═1e4d8b70-c322-4d79-82b7-a6209619ff88
