### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ 321d7ae4-1b91-479d-9c8d-971928d05104
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 0ae54bb6-3f49-11ed-3718-1df4d8077426
using LinearAlgebra,Plots,JLD2,Statistics,PlutoUI

# ╔═╡ 6b0d0352-6f99-42fd-afa4-2abd178c083e
TableOfContents()

# ╔═╡ 5d3b83f3-01b8-49c5-b55e-0067fadce18a
function combine_files(files,path,new_file)
	jobids = Vector{String}(undef,length(files))
	params = Vector{Any}(undef,length(files))
	data = Vector{Array{Float64,3}}(undef,length(files))
	positiondata = Vector{PositionData}(undef, length(files))
	for (i,f) in enumerate(files)
		jobids[i] = load(path*f,"jobid")
		params[i] = load(path*f,"params")
		positiondata[i] = load(path*f,"positiondata")
		data[i] = load(path*f,"data")
	end
	jldopen(path*new_file, "w") do file
        file["data"] = data
		file["positiondata"] = positiondata
        file["params"] = params
        file["jobid"] = jobids
    end
end


# ╔═╡ 7a743bd1-687e-4665-9443-055df6b3df30
path = pwd()*"/data/xx"

# ╔═╡ 553015b5-ae1d-4928-bfa6-1b486a9d2e6e
folders = ["/noisy_chain_pbc_1.9/","/noisy_chain_pbc_1.4/","/noisy_chain_pbc_1.0/","/noisy_chain_pbc_0.8/","/box_pbc_0.8/","/box_pbc_0.5/","/box_pbc_0.1/"]

# ╔═╡ ceea6703-975e-497b-a103-ca80d40b722c
begin
	for (i,folder) in enumerate(folders)
		l = []
		for s in readdir(path*folders[i])
			if !occursin("otoc",s) && !occursin("combined",s)
			append!(l,[s])
			end
		end
		print(l)
		print("\n")
		combine_files(l,path*folder,"combined_file.jld2")
	end
end

# ╔═╡ Cell order:
# ╠═0ae54bb6-3f49-11ed-3718-1df4d8077426
# ╠═321d7ae4-1b91-479d-9c8d-971928d05104
# ╠═6b0d0352-6f99-42fd-afa4-2abd178c083e
# ╠═5d3b83f3-01b8-49c5-b55e-0067fadce18a
# ╠═7a743bd1-687e-4665-9443-055df6b3df30
# ╠═553015b5-ae1d-4928-bfa6-1b486a9d2e6e
# ╠═ceea6703-975e-497b-a103-ca80d40b722c
