### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ de1f56d3-b5b8-4f85-b7fd-3f216e77261d
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ f87c4dd4-9c18-4cf5-9d04-7321841cfe5f
using LinearAlgebra,Plots,JLD2,Statistics,PlutoUI

# ╔═╡ dde1d288-6f3e-47cf-9586-179bae5b5f8a
TableOfContents()

# ╔═╡ 7bdb9e5f-12e1-4277-8bb0-02ce33b41bc4
function combine_files(files,path,new_file)
	jobids = Vector{String}(undef,length(files))
	params = Vector{Any}(undef,length(files))
	data = Vector{Array{Float64,3}}(undef,length(files))
	positiondata = Vector{PositionData}(undef, length(files))
	for (i,f) in enumerate(files)
		jobids[i] = load(path*f,"jobid")
		params[i] = load(path*f,"params")
		positiondata[i] = load(path*f,"positiondata")
		data[i] = load(path*f,"data")[:,:,:]
	end
	jldopen(path*new_file, "w") do file
        file["data"] = data
		file["positiondata"] = positiondata
        file["params"] = params
        file["jobid"] = jobids
    end
end

# ╔═╡ 9257d7cb-a5e6-4e9a-8105-23ccf9d9177c
ρs = [0.5,0.6,0.7,0.8,0.9,1.0,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95]

# ╔═╡ 0cdfc9e9-dd57-4802-ad92-37ccbbb7dcb4
path = pwd() * "/data/"

# ╔═╡ 4e961eea-dacf-431d-aa58-5f9d74c0dfcf
folders = ["noisy_chain_pbc_0.5/","noisy_chain_pbc_0.6/","noisy_chain_pbc_0.7/","noisy_chain_pbc_0.8/","noisy_chain_pbc_0.9/","noisy_chain_pbc_1.00/","noisy_chain_pbc_1.05/","noisy_chain_pbc_1.10/","noisy_chain_pbc_1.15/","noisy_chain_pbc_1.20/","noisy_chain_pbc_1.25/","noisy_chain_pbc_1.30/","noisy_chain_pbc_1.35/","noisy_chain_pbc_1.40/","noisy_chain_pbc_1.45/","noisy_chain_pbc_1.50/","noisy_chain_pbc_1.55/","noisy_chain_pbc_1.60/","noisy_chain_pbc_1.65/","noisy_chain_pbc_1.70/","noisy_chain_pbc_1.75/","noisy_chain_pbc_1.80/","noisy_chain_pbc_1.85/","noisy_chain_pbc_1.90/","noisy_chain_pbc_1.95/"]

# ╔═╡ e0c8f820-25af-472f-9ed0-cf2a78c88af4
begin
	for (i,folder) in enumerate(folders)
		ρ = ρs[i]
		print("---"*folder*"---\n\n")
			l = []
			for s in readdir(path*folders[i])
				if !occursin("otoc",s) && !occursin("combined",s)
					append!(l,[s])
				end
			end
			print("ρ$ρ:\n")
			print(l)
			print("\n\n")
			combine_files(l,path*folder,"combined_file_ρ$ρ.jld2")
		end
end

# ╔═╡ Cell order:
# ╠═f87c4dd4-9c18-4cf5-9d04-7321841cfe5f
# ╠═de1f56d3-b5b8-4f85-b7fd-3f216e77261d
# ╠═dde1d288-6f3e-47cf-9586-179bae5b5f8a
# ╠═7bdb9e5f-12e1-4277-8bb0-02ce33b41bc4
# ╠═9257d7cb-a5e6-4e9a-8105-23ccf9d9177c
# ╠═0cdfc9e9-dd57-4802-ad92-37ccbbb7dcb4
# ╠═4e961eea-dacf-431d-aa58-5f9d74c0dfcf
# ╠═e0c8f820-25af-472f-9ed0-cf2a78c88af4
