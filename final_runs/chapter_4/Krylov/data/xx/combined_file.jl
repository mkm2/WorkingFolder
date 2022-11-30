### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ 11f2d592-3ddb-455b-8774-1982cf079043
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 41a65f30-35c3-4a65-9f80-6ff44d7c07c0
using LinearAlgebra,Plots,JLD2,Statistics,PlutoUI

# ╔═╡ 1bb42fcc-30bf-4fdb-abef-0fec0fdbe4d5
TableOfContents()

# ╔═╡ 3149d476-40e0-4f5f-8e9a-94964e73bd49
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

# ╔═╡ eb17e000-c871-4652-9a57-1ac69f4d7e2a
hs = [3,12,24]

# ╔═╡ 97222b74-f5ae-4e9d-9c94-d3e3f020adf2
path = pwd()

# ╔═╡ acf77651-f95b-4aee-9ff4-4512d57f2120
folders = ["/h = 3/","/h = 12/","/h = 24/"]

# ╔═╡ fc5b4b41-86a7-4db8-b05f-9deafea1038c
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

# ╔═╡ d7ca0087-a858-4ac7-a355-a2899d03a242


# ╔═╡ Cell order:
# ╠═41a65f30-35c3-4a65-9f80-6ff44d7c07c0
# ╠═11f2d592-3ddb-455b-8774-1982cf079043
# ╠═1bb42fcc-30bf-4fdb-abef-0fec0fdbe4d5
# ╠═3149d476-40e0-4f5f-8e9a-94964e73bd49
# ╠═eb17e000-c871-4652-9a57-1ac69f4d7e2a
# ╠═97222b74-f5ae-4e9d-9c94-d3e3f020adf2
# ╠═acf77651-f95b-4aee-9ff4-4512d57f2120
# ╠═fc5b4b41-86a7-4db8-b05f-9deafea1038c
# ╠═d7ca0087-a858-4ac7-a355-a2899d03a242
