### A Pluto.jl notebook ###
# v0.19.15

using Markdown
using InteractiveUtils

# ╔═╡ 87c6500a-5c7f-4023-b1f4-9dc80e091959
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 08de38af-b45f-426e-b24e-d9e219c4c6cf
using LinearAlgebra,Plots,JLD2,Statistics,PlutoUI

# ╔═╡ 81ded523-a2aa-49b7-ad60-95959ef1bc25
TableOfContents()

# ╔═╡ 282a29f3-2b77-47a8-b99d-6bc831295420
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

# ╔═╡ daed97e0-d489-4412-933d-3dabb07713fc
hs = [3,6,9,12]

# ╔═╡ 3319d577-02c1-44e6-82f1-a7e43c690148
path = pwd()

# ╔═╡ a5baa852-1864-4ab4-9f11-ce89ad767ca6
folders = ["/h = 3/","/h = 6/","/h = 9/","/h = 12/"]

# ╔═╡ 692d4636-59d1-4ad9-b18f-bbbee35fdab3
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

# ╔═╡ b9d56792-dc10-4c6e-836f-2b4f27fef0d8
load(path*"/h = 3/combined_file_h3.jld2")

# ╔═╡ 290b0546-706a-4245-86f5-9ad0bb4bd48c
load(path*"/h = 3/7433343_N13_ED.jld2")

# ╔═╡ Cell order:
# ╠═08de38af-b45f-426e-b24e-d9e219c4c6cf
# ╠═87c6500a-5c7f-4023-b1f4-9dc80e091959
# ╠═81ded523-a2aa-49b7-ad60-95959ef1bc25
# ╠═282a29f3-2b77-47a8-b99d-6bc831295420
# ╠═daed97e0-d489-4412-933d-3dabb07713fc
# ╠═3319d577-02c1-44e6-82f1-a7e43c690148
# ╠═a5baa852-1864-4ab4-9f11-ce89ad767ca6
# ╠═692d4636-59d1-4ad9-b18f-bbbee35fdab3
# ╠═b9d56792-dc10-4c6e-836f-2b4f27fef0d8
# ╠═290b0546-706a-4245-86f5-9ad0bb4bd48c
