### A Pluto.jl notebook ###
# v0.19.15

using Markdown
using InteractiveUtils

# ╔═╡ 44ba03b4-e321-48eb-973c-e96ca96934b9
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 3ff499b7-2745-4d7c-9a87-75fd5281e631
using LinearAlgebra,Plots,JLD2,Statistics,PlutoUI

# ╔═╡ fc3bd76c-8b0b-48e9-a3ea-a6bc50f11294
TableOfContents()

# ╔═╡ bd94ca1c-caa8-48cf-84e0-cce369bc86be
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

# ╔═╡ e285de5e-6ccc-492f-a202-de12f9e95043
hs = [12]

# ╔═╡ 673f361a-09a5-42e8-a20b-1ecb58b41567
path = pwd()

# ╔═╡ 9175b59d-3c73-402e-87c0-5bd21c87fdd7
folders = ["/h = 12/"]

# ╔═╡ 2a2a7140-d11c-44bd-bdde-40c5ce3fdf78
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
# ╠═3ff499b7-2745-4d7c-9a87-75fd5281e631
# ╠═44ba03b4-e321-48eb-973c-e96ca96934b9
# ╠═fc3bd76c-8b0b-48e9-a3ea-a6bc50f11294
# ╠═bd94ca1c-caa8-48cf-84e0-cce369bc86be
# ╠═e285de5e-6ccc-492f-a202-de12f9e95043
# ╠═673f361a-09a5-42e8-a20b-1ecb58b41567
# ╠═9175b59d-3c73-402e-87c0-5bd21c87fdd7
# ╠═2a2a7140-d11c-44bd-bdde-40c5ce3fdf78
