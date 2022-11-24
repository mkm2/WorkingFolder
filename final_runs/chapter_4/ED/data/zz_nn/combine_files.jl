### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ 86efe950-ecfe-4b48-bd6d-1e14412a8993
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 041b6699-b8c9-4ddf-b62c-53484d8da45e
using LinearAlgebra,Plots,JLD2,Statistics,PlutoUI

# ╔═╡ 54879234-bccd-406c-b132-82619ea4802b
TableOfContents()

# ╔═╡ 64d6e1cf-b6b4-4205-8d00-5192203c8dd1
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

# ╔═╡ 52b87a98-0de5-4770-9487-ccd36dc3bf2a
hs = [0,3,6,9,12,24]

# ╔═╡ 6024803f-43b4-4e9f-861c-ba093f6c5bf4
path = pwd()

# ╔═╡ 48b87947-ce25-4775-9493-c10cb782bb79
folders = ["/h = 0/","/h = 3/","/h = 6/","/h = 9/","/h = 12/","/h = 24/"]

# ╔═╡ 4c9a4768-85fb-4d3d-bfca-e29ee436a5f4
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
# ╠═041b6699-b8c9-4ddf-b62c-53484d8da45e
# ╠═86efe950-ecfe-4b48-bd6d-1e14412a8993
# ╠═54879234-bccd-406c-b132-82619ea4802b
# ╠═64d6e1cf-b6b4-4205-8d00-5192203c8dd1
# ╠═52b87a98-0de5-4770-9487-ccd36dc3bf2a
# ╠═6024803f-43b4-4e9f-861c-ba093f6c5bf4
# ╠═48b87947-ce25-4775-9493-c10cb782bb79
# ╠═4c9a4768-85fb-4d3d-bfca-e29ee436a5f4
