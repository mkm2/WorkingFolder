### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ f01ab21f-5cca-42fb-ac3a-273a2671f327
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ bf57e49e-2fff-472c-a98c-68fae6d48eaa
using LinearAlgebra,Plots,JLD2,Statistics,PlutoUI

# ╔═╡ d00284ce-51e0-44c7-994e-768562609a95
TableOfContents()

# ╔═╡ cb668854-56c6-495b-846b-fc951e7364bf
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

# ╔═╡ 726123b4-8b7e-48be-9fe0-21c897dd5c74
hs = [0,3,6,9,12,24]

# ╔═╡ 747f34f2-0cf3-478d-bd3a-74cd15a7ce61
path = pwd()

# ╔═╡ 34fe50c5-ef1c-4836-bd9f-791cbc083792
folders = ["/h = 0/","/h = 3/","/h = 6/","/h = 9/","/h = 12/","/h = 24/"]

# ╔═╡ 5ae86193-e973-4276-9d73-3c34fec8ee14
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

# ╔═╡ eadede4c-e406-407d-b372-2a2b498a8158


# ╔═╡ ff6875b2-d2e8-43a5-8244-8d272ef53278


# ╔═╡ bd91c812-dd49-4a28-9af0-c23100193ff6


# ╔═╡ f88acded-dc54-4fb2-9b0e-cf4b476bdfe0


# ╔═╡ 0a568f4f-d062-48ce-9007-986fc5cd0836


# ╔═╡ 3d028312-1f25-4366-a40e-23eb9de86926


# ╔═╡ Cell order:
# ╠═bf57e49e-2fff-472c-a98c-68fae6d48eaa
# ╠═f01ab21f-5cca-42fb-ac3a-273a2671f327
# ╠═d00284ce-51e0-44c7-994e-768562609a95
# ╠═cb668854-56c6-495b-846b-fc951e7364bf
# ╠═726123b4-8b7e-48be-9fe0-21c897dd5c74
# ╠═747f34f2-0cf3-478d-bd3a-74cd15a7ce61
# ╠═34fe50c5-ef1c-4836-bd9f-791cbc083792
# ╠═5ae86193-e973-4276-9d73-3c34fec8ee14
# ╠═eadede4c-e406-407d-b372-2a2b498a8158
# ╠═ff6875b2-d2e8-43a5-8244-8d272ef53278
# ╠═bd91c812-dd49-4a28-9af0-c23100193ff6
# ╠═f88acded-dc54-4fb2-9b0e-cf4b476bdfe0
# ╠═0a568f4f-d062-48ce-9007-986fc5cd0836
# ╠═3d028312-1f25-4366-a40e-23eb9de86926
