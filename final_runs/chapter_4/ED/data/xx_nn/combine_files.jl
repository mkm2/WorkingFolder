### A Pluto.jl notebook ###
# v0.19.17

using Markdown
using InteractiveUtils

# ╔═╡ 3a5ecf44-8289-4488-b6e9-85a37cfab5de
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 6c275e97-eb15-4eae-8683-5c81d88120e9
using LinearAlgebra,Plots,JLD2,Statistics,PlutoUI

# ╔═╡ 0fdcdb44-2d82-4bc2-8e8f-b038477c810b
TableOfContents()

# ╔═╡ 99f98fdb-ea1e-4d33-8e97-087a6e92339f
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

# ╔═╡ 1ec9e5f3-9d38-466e-8c92-875e8ed61963
hs = [0,3,6,9,12,24,24,50,75,50,75,100,24]

# ╔═╡ 59d6b2eb-edb0-4238-b1d6-0358feca23ee
path = pwd()

# ╔═╡ 027f13d1-368d-4233-b84d-f8cc1c487910
folders = ["/h = 0/","/h = 3/","/h = 6/","/h = 9/","/h = 12/","/h = 24/","/h = 24 ext/","/h = 50/","/h = 75/","/h = 50 long/","/h = 75 long/","/h = 100 long/","/h = 24 long/"]

# ╔═╡ d4ce4c4d-d791-4c5a-b30c-3382e0c780c2
begin
	for (i,folder) in enumerate(folders)
		h = hs[i]
		print("---"*folder*"---\n\n")
			l = []
			for s in readdir(path*folders[i])
				if !occursin("otoc",s) && !occursin("combined",s)
					append!(l,[s])
				end
			end
			print("h$h:\n")
			print(l)
			print("\n\n")
			combine_files(l,path*folder,"combined_file_h$h.jld2")
		end
end

# ╔═╡ Cell order:
# ╠═6c275e97-eb15-4eae-8683-5c81d88120e9
# ╠═3a5ecf44-8289-4488-b6e9-85a37cfab5de
# ╠═0fdcdb44-2d82-4bc2-8e8f-b038477c810b
# ╠═99f98fdb-ea1e-4d33-8e97-087a6e92339f
# ╠═1ec9e5f3-9d38-466e-8c92-875e8ed61963
# ╠═59d6b2eb-edb0-4238-b1d6-0358feca23ee
# ╠═027f13d1-368d-4233-b84d-f8cc1c487910
# ╠═d4ce4c4d-d791-4c5a-b30c-3382e0c780c2
