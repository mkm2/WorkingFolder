### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ f76a6752-ccd3-4601-8498-916de97c8a4f
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 300894b6-3cf1-11ed-23f6-8f618c96822a
using LinearAlgebra,Plots,JLD2,Statistics,PlutoUI

# ╔═╡ 729d35fd-db62-4f26-8b7b-6aea5733b940
using SparseArrays

# ╔═╡ 6267285a-85f0-4c8c-859a-f0e9ae2b2dcf
TableOfContents()

# ╔═╡ 649d4204-cb0f-4126-996b-422286123841
function combine_files(files,path,new_file)
	jobids = Vector{String}(undef,length(files))
	params = Vector{Any}(undef,length(files))
	data = Vector{Array{Float64,3}}(undef,length(files))
	for (i,f) in enumerate(files)
		jobids[i] = load(path*f,"jobid")
		params[i] = load(path*f,"params")
		data[i] = load(path*f,"data")
	end
	jldopen(path*new_file, "w") do file
        file["data"] = data
        file["params"] = params
        file["jobid"] = jobids
    end
end

# ╔═╡ 0719d385-9603-45f7-8da3-93a1ac29f8f5
begin
		function disorder_mean(A,n_shots)
			return mean(A[:,:,1:n_shots];dims=3)[:,:,1]
		end
		
		function disorder_std(A,n_shots)
			return std(A[:,:,1:n_shots];dims=3)[:,:,1]
		end
	
		#function pos_mean(A)
		#	return mean(A;dims=2)[:,1]
		#end
		
		#function reduce_by_last(A)
			#return A[:,:,1]
		#end
end

# ╔═╡ 4b5386b8-79bf-46b2-943f-6463f2f7b97e
trange = logrange(-2,10,1e10)

# ╔═╡ d0a26952-c154-4aa1-b963-982ac5b50d4d
length(trange)

# ╔═╡ 3cd6b1e3-a2d5-45c2-867c-e44c2233538e
md"# tmax = 1e10"

# ╔═╡ c8c6fd98-f5c9-4bb4-86d6-af8a7456eaa6
md"## h = 0"

# ╔═╡ 44795427-654c-43cb-9c7b-99e56d8f4efa
N = 15

# ╔═╡ cec3b3a8-fd64-4922-9705-cf0e8c76c13c
path = pwd()*"/"

# ╔═╡ 8dd6f239-f1be-4d15-b967-f9df3c02f820


# ╔═╡ afe723da-0e2a-4a2b-9684-a8e3919fa29c


# ╔═╡ d6a2c47b-2bc2-486c-a4ff-2676e8488cb1
md"## h = 4"

# ╔═╡ 42acf62e-a8fe-44cf-a895-0b3fb5274aec
shots4 = 100

# ╔═╡ 1017c624-4279-4b96-bc3f-90a5bfe16576
begin
	f_4 = "7143487-7143488_N15_ED.jld2"
	combine_files(["7143487_N15_ED.jld2","7143488_N15_ED.jld2"],path,f_4)
	
	jobids_4 = load(path*f_4,"jobid")
	params_4 = load(path*f_4,"params")
	data_4 = 2*ones(110,N,100)-2*cat(load(path*f_4,"data")...,dims=3)

	size(data_4)
end

# ╔═╡ 6bc11a97-315b-49df-990d-718701e717f1
begin
	data_4_mean = disorder_mean(data_4,size(data_4)[3])
	data_4_std = disorder_std(data_4,size(data_4)[3])
end

# ╔═╡ bbb0bd91-e295-4e81-b38c-a37bcf00f3fe
plot(trange[2:110],data_4_mean[2:110,:],xaxis=:log,legend=nothing,ribbon=data_4_std[2:110,:]/sqrt(shots4))

# ╔═╡ 8207f588-ce60-4165-aa1e-7efb26937252
heatmap(1:15,trange[2:110],data_4_mean[2:110,:],yaxis=:log)

# ╔═╡ 2cc9f257-fc72-4a10-a543-f088d7a925ae


# ╔═╡ 239c9d65-5a98-4899-9cac-61b09e7eaa3e


# ╔═╡ 058d0476-cb06-44ab-9f9a-bbc92a71b096
md"## h = 12"

# ╔═╡ b9dea71c-af5b-4e5f-b925-2b268f1944b0
shots12 = 100

# ╔═╡ ea17c8f7-59c6-4fa0-b205-9a01272f2a69
begin
	f_12 = "7143489-7143490_N15_ED.jld2"
	combine_files(["7143489_N15_ED.jld2","7143490_N15_ED.jld2"],path,f_12)
	
	jobids_12 = load(path*f_12,"jobid")
	params_12 = load(path*f_12,"params")
	data_12 = 2*ones(110,N,100)-2*cat(load(path*f_12,"data")...,dims=3)

	size(data_12)
end

# ╔═╡ 65c3da40-649f-41ce-aabd-075136ab2c66
begin
	data_12_mean = disorder_mean(data_12,shots12)
	data_12_std = disorder_std(data_12,shots12)
end

# ╔═╡ 88046dbf-0921-484d-88a5-0e13abae4e9f
plot(trange[2:110],data_12_mean[2:110,:],xaxis=:log,legend=nothing,ribbon=data_12_std[2:110,:]/sqrt(shots12))

# ╔═╡ b70c4906-fc45-4a67-98fb-d92d6b7c9b0b
heatmap(1:15,trange[2:110],data_12_mean[2:110,:],yaxis=:log)

# ╔═╡ f041788a-d9c9-49dc-a9f9-1e2404fb866f
md"## h = 20"

# ╔═╡ 7584be28-2fc1-4de4-ab94-feaa6d8f62e2
shots20 = 100

# ╔═╡ 99488b4f-7d12-4cf7-84a2-c2d7e20794f4
begin
	f_20 = "7143491-7143492_N15_ED.jld2"
	combine_files(["7143491_N15_ED.jld2","7143492_N15_ED.jld2"],path,f_20)
	
	jobids_20 = load(path*f_20,"jobid")
	params_20 = load(path*f_20,"params")
	data_20 = 2*ones(110,N,100)-2*cat(load(path*f_20,"data")...,dims=3)

	size(data_20)
end

# ╔═╡ fe85d320-14c8-4c63-8431-a45439c0141d
begin
	data_20_mean = disorder_mean(data_20,shots20)
	data_20_std = disorder_std(data_20,shots20)
end

# ╔═╡ dfeb47d9-b990-4b72-8f3d-6768c562f547
plot(trange[2:110],data_20_mean[2:110,:],xaxis=:log,legend=nothing,ribbon=data_20_std[2:110,:]/sqrt(shots20))

# ╔═╡ 3c6f6533-9a2c-4663-a6c3-0fe813ac9a0c
heatmap(1:15,trange[2:110],data_20_mean[2:110,:],yaxis=:log)

# ╔═╡ 140693d2-1d1e-4445-bd0a-3aaa1fde6ce4


# ╔═╡ 46a58d6f-90ba-45a8-b03e-4a0aff3ebaa1
M = 8

# ╔═╡ 851d9d59-a662-4946-9ad1-2cca54ecb234
H = xxz(M)

# ╔═╡ fc838549-7ced-4ffb-9c4d-797ad8398eed
A = convert(SparseMatrixCSC{ComplexF64,Int64},single_spin_op(σz,2,M))

# ╔═╡ af4eb76f-2933-4fcb-aa8b-caba4e29b9ec
begin
	B = σz
	B = convert(SparseMatrixCSC{ComplexF64,Int64},B)
end

# ╔═╡ a86c6b6d-19c0-4663-a9c6-55ad7dca90d0
s = 10

# ╔═╡ a98acd99-7990-4b76-9b33-df31f652fb32
ψ = convert(Vector{ComplexF64},neel_state(M,2^M))

# ╔═╡ 2c5e0641-8e8a-4ea4-ae71-5076bd93601d


# ╔═╡ d5528466-001e-4e16-a320-fd5007879959
test = Diag_OTOCψ(Matrix(H),A,B,trange,ψ)

# ╔═╡ 083adbec-7e1f-4bff-b3c4-8a37ff7c8de3
plot(trange[2:110],2*(ones(109,11)-test[2:110,:]),xaxis=:log,legend=nothing)

# ╔═╡ e445e317-eff1-40e1-8cec-26d5337755ad


# ╔═╡ Cell order:
# ╠═300894b6-3cf1-11ed-23f6-8f618c96822a
# ╠═f76a6752-ccd3-4601-8498-916de97c8a4f
# ╠═6267285a-85f0-4c8c-859a-f0e9ae2b2dcf
# ╠═649d4204-cb0f-4126-996b-422286123841
# ╠═0719d385-9603-45f7-8da3-93a1ac29f8f5
# ╠═4b5386b8-79bf-46b2-943f-6463f2f7b97e
# ╠═d0a26952-c154-4aa1-b963-982ac5b50d4d
# ╠═3cd6b1e3-a2d5-45c2-867c-e44c2233538e
# ╠═c8c6fd98-f5c9-4bb4-86d6-af8a7456eaa6
# ╠═44795427-654c-43cb-9c7b-99e56d8f4efa
# ╠═cec3b3a8-fd64-4922-9705-cf0e8c76c13c
# ╠═8dd6f239-f1be-4d15-b967-f9df3c02f820
# ╠═afe723da-0e2a-4a2b-9684-a8e3919fa29c
# ╠═d6a2c47b-2bc2-486c-a4ff-2676e8488cb1
# ╠═42acf62e-a8fe-44cf-a895-0b3fb5274aec
# ╠═1017c624-4279-4b96-bc3f-90a5bfe16576
# ╠═6bc11a97-315b-49df-990d-718701e717f1
# ╠═bbb0bd91-e295-4e81-b38c-a37bcf00f3fe
# ╠═8207f588-ce60-4165-aa1e-7efb26937252
# ╠═2cc9f257-fc72-4a10-a543-f088d7a925ae
# ╠═239c9d65-5a98-4899-9cac-61b09e7eaa3e
# ╠═058d0476-cb06-44ab-9f9a-bbc92a71b096
# ╠═b9dea71c-af5b-4e5f-b925-2b268f1944b0
# ╠═ea17c8f7-59c6-4fa0-b205-9a01272f2a69
# ╠═65c3da40-649f-41ce-aabd-075136ab2c66
# ╠═88046dbf-0921-484d-88a5-0e13abae4e9f
# ╠═b70c4906-fc45-4a67-98fb-d92d6b7c9b0b
# ╠═f041788a-d9c9-49dc-a9f9-1e2404fb866f
# ╠═7584be28-2fc1-4de4-ab94-feaa6d8f62e2
# ╠═99488b4f-7d12-4cf7-84a2-c2d7e20794f4
# ╠═fe85d320-14c8-4c63-8431-a45439c0141d
# ╠═dfeb47d9-b990-4b72-8f3d-6768c562f547
# ╠═3c6f6533-9a2c-4663-a6c3-0fe813ac9a0c
# ╠═140693d2-1d1e-4445-bd0a-3aaa1fde6ce4
# ╠═46a58d6f-90ba-45a8-b03e-4a0aff3ebaa1
# ╠═851d9d59-a662-4946-9ad1-2cca54ecb234
# ╠═729d35fd-db62-4f26-8b7b-6aea5733b940
# ╠═fc838549-7ced-4ffb-9c4d-797ad8398eed
# ╠═af4eb76f-2933-4fcb-aa8b-caba4e29b9ec
# ╠═a86c6b6d-19c0-4663-a9c6-55ad7dca90d0
# ╠═a98acd99-7990-4b76-9b33-df31f652fb32
# ╠═2c5e0641-8e8a-4ea4-ae71-5076bd93601d
# ╠═d5528466-001e-4e16-a320-fd5007879959
# ╠═083adbec-7e1f-4bff-b3c4-8a37ff7c8de3
# ╠═e445e317-eff1-40e1-8cec-26d5337755ad
