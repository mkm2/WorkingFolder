### A Pluto.jl notebook ###
# v0.19.15

using Markdown
using InteractiveUtils

# ╔═╡ 8e640669-b5f7-4035-a7a1-d3ef9bed06d2
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 0fd6acd6-f8e5-4b3e-a657-0b17dc8ec4e9
using LinearAlgebra,Plots,JLD2,Statistics,PlutoUI

# ╔═╡ 28f55a00-3ae8-4475-8e1f-5e4b653a0e81
TableOfContents()

# ╔═╡ 079227bb-fb02-4df5-90fd-cba811ea7a34
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

# ╔═╡ f3c6d15e-ad24-4d43-aa0f-efc10b9eb76b
ds = [1.99,1.8,1.5,1.2,0.9,0.6]

# ╔═╡ 2a5697cc-83e6-4687-8047-0db6cb6ba2d2
path = pwd()

# ╔═╡ 1bb98880-6a73-47fe-88eb-d00afe929e0d
folders = ["/d=1.99/","/d=1.8/","/d=1.5/","/d=1.2/","/d=0.9/","/d=0.6/"]

# ╔═╡ 5c11434b-3971-4b06-92ab-247559ee787b
begin
	for (i,folder) in enumerate(folders)
		d = ds[i]
		print("---"*folder*"---\n\n")
			l = []
			for s in readdir(path*folders[i])
				if !occursin("otoc",s) && !occursin("combined",s)
					append!(l,[s])
				end
			end
			print("d$d:\n")
			print(l)
			print("\n\n")
			combine_files(l,path*folder,"combined_file_d$d.jld2")
		end
end

# ╔═╡ Cell order:
# ╠═0fd6acd6-f8e5-4b3e-a657-0b17dc8ec4e9
# ╠═8e640669-b5f7-4035-a7a1-d3ef9bed06d2
# ╠═28f55a00-3ae8-4475-8e1f-5e4b653a0e81
# ╠═079227bb-fb02-4df5-90fd-cba811ea7a34
# ╠═f3c6d15e-ad24-4d43-aa0f-efc10b9eb76b
# ╠═2a5697cc-83e6-4687-8047-0db6cb6ba2d2
# ╠═1bb98880-6a73-47fe-88eb-d00afe929e0d
# ╠═5c11434b-3971-4b06-92ab-247559ee787b
