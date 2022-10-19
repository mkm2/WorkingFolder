### A Pluto.jl notebook ###
# v0.19.13

using Markdown
using InteractiveUtils

# ╔═╡ b2e67622-89ae-4292-9937-a7c885d7059c
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 92b1c6ac-d4ec-4ab6-a39b-b94a04b45443
using LinearAlgebra,Plots,JLD2,Statistics,PlutoUI

# ╔═╡ 65ef8145-97cb-4c67-be3a-ceb5d4cba84e
TableOfContents()

# ╔═╡ d82b5f24-3765-431d-9def-d7d7d451df19
function combine_files(files,path,new_file)
	jobids = Vector{String}(undef,length(files))
	params = Vector{Any}(undef,length(files))
	data = Vector{Array{Float64,4}}(undef,length(files))
	#positiondata = Vector{PositionData}(undef, length(files))
	for (i,f) in enumerate(files)
		jobids[i] = load(path*f,"jobid")
		params[i] = load(path*f,"params")
		#positiondata[i] = load(path*f,"positiondata")
		data[i] = load(path*f,"data")[1:31,:,:,:]
	end
	jldopen(path*new_file, "w") do file
        file["data"] = data
		#file["positiondata"] = positiondata
        file["params"] = params
        file["jobid"] = jobids
    end
end

# ╔═╡ 39e12935-fea0-462b-bbfc-cc0c79f8a0ea
Ns = [13,15,17,19,21]

# ╔═╡ 2cf67ce8-f61e-451c-a02a-e959adea2b3c


# ╔═╡ 14ce3d7c-6bcd-4444-be15-b63349911134
string(13)

# ╔═╡ 0b436632-b156-417d-8f69-ea9c8d228d98
path = pwd()

# ╔═╡ bd284fd3-5117-40e9-b2e2-2d02e350d3c5
folders = ["/BS/","/RPS/","/RS/"]

# ╔═╡ c92d108e-06b3-48ec-b542-404d4e432e04
begin
	for (i,folder) in enumerate(folders)
		print("---"*folder*"---\n\n")
		for N in Ns
			l = []
			for s in readdir(path*folders[i])
				if !occursin("otoc",s) && !occursin("combined",s) && occursin("N"*string(N),s)
					append!(l,[s])
				end
			end
			print("N$N:\n")
			print(l)
			print("\n\n")
			combine_files(l,path*folder,"combined_file_N$N.jld2")
		end
		
	end
end

# ╔═╡ Cell order:
# ╠═92b1c6ac-d4ec-4ab6-a39b-b94a04b45443
# ╠═b2e67622-89ae-4292-9937-a7c885d7059c
# ╠═65ef8145-97cb-4c67-be3a-ceb5d4cba84e
# ╠═d82b5f24-3765-431d-9def-d7d7d451df19
# ╠═39e12935-fea0-462b-bbfc-cc0c79f8a0ea
# ╠═2cf67ce8-f61e-451c-a02a-e959adea2b3c
# ╠═14ce3d7c-6bcd-4444-be15-b63349911134
# ╠═0b436632-b156-417d-8f69-ea9c8d228d98
# ╠═bd284fd3-5117-40e9-b2e2-2d02e350d3c5
# ╠═c92d108e-06b3-48ec-b542-404d4e432e04
