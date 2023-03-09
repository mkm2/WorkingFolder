### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ 51a5517a-cadb-4112-af62-019a3291a575
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 8cda2560-b2e9-11ed-3cff-f588418f9ea2
using LinearAlgebra,Plots,JLD2,Statistics,PlutoUI

# ╔═╡ d841a4c6-2bb0-493e-a994-043203f0e11f
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

# ╔═╡ 825b33b7-97e2-4283-9979-c118615c82bb
path = pwd()

# ╔═╡ c291e688-5ba9-4be2-a081-1e8e20a7ce68
hs = [0,3,6,9,12,15,18,21,24]

# ╔═╡ 5f906c64-65fd-42cd-ad2e-35fe0b2c4aa5
begin
	top_folders = ["/alpha = 3/PBC","/alpha = 3/OBC","/alpha = 6/PBC","/alpha = 6/OBC","/nn/PBC","/nn/OBC"]
	h_folders = ["/h = 0/","/h = 3/","/h = 6/","/h = 9/","/h = 12/","/h = 15/","/h = 18/","/h = 21/","/h = 24/"]
end

# ╔═╡ 46d584b5-aa73-45f9-80b7-d996a34bde90
begin
	for tf in top_folders
		for (i,hfolder) in enumerate(h_folders)
			h = hs[i]
			print("---"*tf*hfolder*"---\n\n")
				l = []
				for s in readdir(path*tf*hfolder)
					if !occursin("otoc",s) && !occursin("combined",s)
						append!(l,[s])
					end
				end
				print("h$h:\n")
				print(l)
				print("\n\n")
				combine_files(l,path*tf*hfolder,"combined_file_h$h.jld2")
			end
	end
end

# ╔═╡ baa4d46b-385f-4bf5-b14f-04aede079713


# ╔═╡ Cell order:
# ╠═8cda2560-b2e9-11ed-3cff-f588418f9ea2
# ╠═51a5517a-cadb-4112-af62-019a3291a575
# ╠═d841a4c6-2bb0-493e-a994-043203f0e11f
# ╠═825b33b7-97e2-4283-9979-c118615c82bb
# ╠═c291e688-5ba9-4be2-a081-1e8e20a7ce68
# ╠═5f906c64-65fd-42cd-ad2e-35fe0b2c4aa5
# ╠═46d584b5-aa73-45f9-80b7-d996a34bde90
# ╠═baa4d46b-385f-4bf5-b14f-04aede079713
