### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ 1d3bebc6-d0fd-4ae7-b0e7-9451a305c28e
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 3370f337-edb9-4f45-9fb3-9aae676e9ae5
using LinearAlgebra,Plots,JLD2,Statistics,PlutoUI

# ╔═╡ 02efb613-941f-4f29-a098-afed14adb8b6
TableOfContents()

# ╔═╡ 5fac7461-6271-4fe1-bc5f-f014ec3a02c7
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

# ╔═╡ 0f168572-6e04-4cca-8508-53a8bc7ef4a0
names = ["xx","xxnn"]

# ╔═╡ f0738b9f-4d2f-4ec8-a11e-657a32d0ec12
path = pwd()

# ╔═╡ b9b2d082-b805-4606-ae52-4f4b72e3a1df
folders = ["/OBC_offcenter/xx_h24/","/OBC_offcenter/xx_nn_h24/"]

# ╔═╡ a83cd9d9-470d-464c-a909-900925e68933
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
# ╠═3370f337-edb9-4f45-9fb3-9aae676e9ae5
# ╠═1d3bebc6-d0fd-4ae7-b0e7-9451a305c28e
# ╠═02efb613-941f-4f29-a098-afed14adb8b6
# ╠═5fac7461-6271-4fe1-bc5f-f014ec3a02c7
# ╠═0f168572-6e04-4cca-8508-53a8bc7ef4a0
# ╠═f0738b9f-4d2f-4ec8-a11e-657a32d0ec12
# ╠═b9b2d082-b805-4606-ae52-4f4b72e3a1df
# ╠═a83cd9d9-470d-464c-a909-900925e68933
