### A Pluto.jl notebook ###
# v0.19.13

using Markdown
using InteractiveUtils

# ╔═╡ db8bcfe9-3519-43c4-abb8-df5a3a757c73
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ b630a7e4-3e73-11ed-1baf-c976bb434b06
using LinearAlgebra,Plots,JLD2,Statistics,PlutoUI

# ╔═╡ 7e135432-4ad6-47c0-a2d6-52da562f66e8
TableOfContents()

# ╔═╡ 166576ac-961d-4a22-9406-327e7e76899d
length(logrange(-2,10,1e10))

# ╔═╡ edb5ce15-8d34-4ff3-ab4a-873dbf1a61a3
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

# ╔═╡ c6561277-446a-49f4-ad4a-a9a60ffe57df
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

# ╔═╡ 3b511bc9-840b-47a7-b4a2-363ca0188e0e
trange = logrange(-2,10,1e10)

# ╔═╡ ecb93b7f-e7bf-4cfd-b0dd-eab8ecb699fe
T = length(trange)

# ╔═╡ 6b87691e-c010-41c2-ac8d-721c6ce972bb
md"# σxσx: N = 13, tmax = 1e10"

# ╔═╡ 71b9bf79-5135-4ee5-bcdc-7ba1642910af
begin
	N = 13
	path = pwd() *"/"
	shots = 100
end

# ╔═╡ 6d7f5c78-8082-491b-a884-5839185641fc
md"# PowerLaw OBC"

# ╔═╡ 6c0d6752-667b-49a7-8f7f-7a039e7104a6
md"## h = 0"

# ╔═╡ 4af95c99-f7a3-4d71-b050-a37e8a5f1f03
begin
	f_0 = "XXX_N13_ED.jld2"
	
	jobids_0 = load(path*f_0,"jobid")
	params_0 = load(path*f_0,"params")
	data_0 = 2*ones(T,N,50)-2*load(path*f_0,"data")

	size(data_0)
end

# ╔═╡ bc844145-f979-4277-a54c-7e902e47ccf0
begin
	data_0_mean = disorder_mean(data_0,shots)
	data_0_std = disorder_std(data_0,shots)
end

# ╔═╡ 4c422f78-8269-491d-9ba3-2df6bfef198d
plot(trange[2:T],data_0_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_0_std[2:T,:]/sqrt(shots))

# ╔═╡ e22ce722-a048-4df1-b77e-97a113defa53
heatmap(1:N,trange[2:T],data_0_mean[2:T,:],yaxis=:log)

# ╔═╡ 9900588f-d02f-4e46-b136-9903da1c121a
md"## h = 4"

# ╔═╡ 88b944a6-ae83-419d-a3b5-988c512bad07
begin
	f_4 = "7180016-7180017_N13_ED.jld2"
	combine_files(["7180016_N13_ED.jld2","7180017_N13_ED.jld2"],path,f_4)
	
	jobids_4 = load(path*f_4,"jobid")
	params_4 = load(path*f_4,"params")
	data_4 = 2*ones(T,N,shots)-2*cat(load(path*f_4,"data")...,dims=3)

	params_4
end

# ╔═╡ fa96cb73-d1cb-401c-addb-d86b8c73dd9a
begin
	data_4_mean = disorder_mean(data_4,shots)
	data_4_std = disorder_std(data_4,shots)
end

# ╔═╡ 9f5dcb52-8418-417b-b7a3-0f7f7682fcbc
plot(trange[2:T],data_4_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_4_std[2:T,:]/sqrt(shots),ylabel="OTOC xx",xlabel="Jt",title="Power Law OBC, h = 4")

# ╔═╡ a148f020-cb39-4820-b868-503d2dca252c
heatmap(1:N,trange[2:T],data_4_mean[2:T,:],yaxis=:log)

# ╔═╡ 10ca2e1c-d5dc-4c39-bb86-0981d951d283
md"## h = 12"

# ╔═╡ 880c16dd-240b-4ad5-92f1-54b923ef75da
begin
	f_12 = "7180018-7180019_N13_ED.jld2"
	combine_files(["7180018_N13_ED.jld2","7180019_N13_ED.jld2"],path,f_12)
	
	jobids_12 = load(path*f_12,"jobid")
	params_12 = load(path*f_12,"params")
	data_12 = 2*ones(T,N,shots)-2*cat(load(path*f_12,"data")...,dims=3)

	params_12
end

# ╔═╡ 8bbe06ac-a2da-4ceb-bc5c-3c442a78c54c
begin
	data_12_mean = disorder_mean(data_12,shots)
	data_12_std = disorder_std(data_12,shots)
end

# ╔═╡ be7d456d-b84d-47fd-a750-949227cbea71
10. .^LinRange(-3,8,120)

# ╔═╡ e8d0bd88-11aa-419a-98bc-05c3cfc3cb13
data_12_mean

# ╔═╡ c4f385eb-6e67-402d-9f98-c8ac8621c9a2
begin
	plot(trange[2:T],data_12_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_12_std[2:T,:]/sqrt(shots),ylabel="OTOC xx",xlabel="Jt",title="Power Law OBC, h= 12",ylims=[0.001,3.2],xlim=[1e-3,1e6])#,xlim=[1e-2,1e8],yaxis=:log)
	#plot!(10. .^ LinRange(-3,8,100),ones(100),xaxis=:log,markershape=:cross)
end

# ╔═╡ 4302e49d-d06a-4e22-86a1-75669096779f
heatmap(1:N,trange[2:T],data_12_mean[2:T,:],yaxis=:log)

# ╔═╡ 9380deb5-848b-40f9-a58d-0a8d157c90b1
md"# NearestNeighbour OBC"

# ╔═╡ 6602f0c5-1445-417a-b2ca-490622c03169
md"## h = 0"

# ╔═╡ bcc17df1-370a-4151-9433-688f3195bc8e
begin
	f_0_nn = "7180031_N13_ED.jld2"
	
	jobids_0_nn = load(path*f_0_nn,"jobid")
	params_0_nn = load(path*f_0_nn,"params")
	data_0_nn = 2*ones(T,N,50)-2*load(path*f_0_nn,"data")

	params_0_nn
end

# ╔═╡ 6f8b67b3-796f-4246-a250-0a023f517ace
begin
	data_0_nn_mean = disorder_mean(data_0_nn,50)
	data_0_nn_std = disorder_std(data_0_nn,50)
end

# ╔═╡ a939797d-8e74-4d0f-bfff-deec14b7ac51
plot(trange[2:T],abs.(data_0_nn_mean[2:110,:]),xaxis=:log,legend=:bottomright,ylabel="OTOC xx",xlabel="Jt",title="Nearest Neighbour OBC, h = 0",yaxis=:linear)

# ╔═╡ 33282538-dadb-4487-986e-ea89f33ebe88
heatmap(1:N,trange[2:T],data_0_nn_mean[2:T,:],yaxis=:log)

# ╔═╡ ac5908dc-61e9-479b-a33c-50e70abd829f
md"## h = 4"

# ╔═╡ 72391c41-6e9a-4c4c-a766-57b27d10445e
begin
	f_4_nn = "7180032-7180033_N13_ED.jld2"
	combine_files(["7180032_N13_ED.jld2","7180033_N13_ED.jld2"],path,f_4_nn)
	
	jobids_4_nn = load(path*f_4_nn,"jobid")
	params_4_nn = load(path*f_4_nn,"params")
	data_4_nn = 2*ones(T,N,shots)-2*cat(load(path*f_4_nn,"data")...,dims=3)

	params_4_nn
end

# ╔═╡ 8679c844-ceb7-44c2-9e88-e13456c77933
begin
	data_4_nn_mean = disorder_mean(data_4_nn,shots)
	data_4_nn_std = disorder_std(data_4_nn,shots)
end

# ╔═╡ 0d72b7be-8a28-400d-9c55-a9ada9d62322
plot(trange[2:T],data_4_nn_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_4_nn_std[2:T,:]/sqrt(shots),ylabel="OTOC xx",xlabel="Jt",title="Nearest Neighbour OBC, h = 4")

# ╔═╡ e7af56ef-24ed-4e49-917e-b05fa0473c0c
heatmap(1:N,trange[2:T],data_4_nn_mean[2:T,:],yaxis=:log)

# ╔═╡ a573b238-1a93-46d1-b7f5-3deccc02bdc5
md"## h = 12"

# ╔═╡ b72312cb-98a8-4a93-bd9f-8a5b5e8c6750
begin
	f_12_nn = "7180034-7180035_N13_ED.jld2"
	combine_files(["7180034_N13_ED.jld2","7180035_N13_ED.jld2"],path,f_12_nn)
	
	jobids_12_nn = load(path*f_12_nn,"jobid")
	params_12_nn = load(path*f_12_nn,"params")
	data_12_nn = 2*ones(T,N,shots)-2*cat(load(path*f_12_nn,"data")...,dims=3)

	params_12_nn
end

# ╔═╡ 88300d49-ecb9-4f08-9262-f109b8297b03
begin
	data_12_nn_mean = disorder_mean(data_12_nn,shots)
	data_12_nn_std = disorder_std(data_12_nn,shots)
end

# ╔═╡ 44181015-34a4-47ef-aefb-6f8509fdfca1
begin
	plot(trange[2:T],abs.(data_12_nn_mean[2:110,:]),xaxis=:log,legend=nothing,ribbon=data_12_nn_std[2:T,:]/sqrt(shots),ylabel="OTOC xx",xlabel="Jt",title="Nearest Neighbour OBC, h = 12",xlim=[1e-2,1e1],yaxis=:log)
	plot!(10. .^ LinRange(-3,1,100),ones(100),xaxis=:log,markershape=:cross)
end

# ╔═╡ 50f76e36-ec8b-4989-bfe2-d109f811ae26
log10(5)

# ╔═╡ 54d6b8bc-6b9b-4243-9aea-56eafbf9f701
heatmap(1:N,trange[2:T],data_12_nn_mean[2:T,:],yaxis=:log)

# ╔═╡ b969969c-967e-4ff6-b781-62087e167f29


# ╔═╡ d21541e3-ba62-4dd5-844e-1751cf134d29
begin
	plot(trange[2:T],data_12_nn_mean[2:110,1:3]-data_12_mean[2:110,1:3],xaxis=:log,legend=:bottomright,ribbon=sqrt.(data_12_nn_std[2:T,:].^2+data_12_std[2:T,:].^2)/sqrt(shots))
end

# ╔═╡ f55a5ebe-b8cf-4772-8aec-ee201db139ef
md"# PowerLaw PBC"

# ╔═╡ 7864764e-8f34-43b9-bedd-970847d38688
md"## h = 0"

# ╔═╡ c9d599c4-6cd8-4733-9603-8bfffd3ac836
trange

# ╔═╡ 9a75ae46-c32d-4d0d-a7d6-e1df83f17129
begin
	f_0_pbc = "7194612_N13_ED.jld2"
	
	jobids_0_pbc = load(path*f_0_pbc,"jobid")
	params_0_pbc = load(path*f_0_pbc,"params")
	data_0_pbc = 2*ones(T,N,50)-2*load(path*f_0_pbc,"data")

	params_0_pbc
end

# ╔═╡ 9dcfe111-e17f-4b7c-807b-a5692d88083a
begin
	data_0_pbc_mean = disorder_mean(data_0_pbc,50)
	data_0_pbc_std = disorder_std(data_0_pbc,50)
end

# ╔═╡ 0a04bf3f-9783-4d1d-99e3-a00e7d844d07
plot(trange[2:T],data_0_pbc_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_0_pbc_std[2:T,:]/sqrt(shots),ylabel="OTOC xx",xlabel="Jt",title="Power Law PBC, h = 0")

# ╔═╡ 5f376697-c573-4871-b8d1-758f9b46b92b
heatmap(1:N,trange[2:T],data_0_pbc_mean[2:T,:],yaxis=:log)

# ╔═╡ 386d9b58-bf17-49d9-a2b8-443aad198c6c
md"## h = 4"

# ╔═╡ 675bca7e-d84a-4774-89f8-9f67ea3e4d2c
begin
	f_4_pbc = "7180021-7180022_N13_ED.jld2"
	combine_files(["7180021_N13_ED.jld2","7180022_N13_ED.jld2"],path,f_4_pbc)
	
	jobids_4_pbc = load(path*f_4_pbc,"jobid")
	params_4_pbc = load(path*f_4_pbc,"params")
	data_4_pbc = 2*ones(T,N,shots)-2*cat(load(path*f_4_pbc,"data")...,dims=3)

	params_4_pbc
end

# ╔═╡ cbbe5b27-2a57-4296-a4ef-b87c0502af4d
begin
	data_4_pbc_mean = disorder_mean(data_4_pbc,shots)
	data_4_pbc_std = disorder_std(data_4_pbc,shots)
end

# ╔═╡ da1670a7-070c-4533-bef5-60b07a3f5edf
plot(trange[2:T],data_4_pbc_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_4_pbc_std[2:T,:]/sqrt(shots),ylabel="OTOC xx",xlabel="Jt",title="Power Law PBC, h= 4")

# ╔═╡ d7e301c3-4e11-444e-80ac-8691b9d4c556
heatmap(1:N,trange[2:T],data_4_pbc_mean[2:T,:],yaxis=:log)

# ╔═╡ 727e60b7-ec48-4a13-b4c7-48b68f27f9e6
md"## h = 12"

# ╔═╡ 96da9674-935c-444a-81ae-a98fa54a95d4
begin
	f_12_pbc = "7194613-7194614_N13_ED.jld2"
	combine_files(["7194613_N13_ED.jld2","7194614_N13_ED.jld2"],path,f_12_pbc)
	
	jobids_12_pbc = load(path*f_12_pbc,"jobid")
	params_12_pbc = load(path*f_12_pbc,"params")
	data_12_pbc = 2*ones(T,N,shots)-2*cat(load(path*f_12_pbc,"data")...,dims=3)

	params_12_pbc
end

# ╔═╡ b2af0bf9-f305-48c6-8295-c560bb9d25ec
begin
	data_12_pbc_mean = disorder_mean(data_12_pbc,shots)
	data_12_pbc_std = disorder_std(data_12_pbc,shots)
end

# ╔═╡ 9378d2bd-2e31-4a7d-a81c-574e70c9a171
plot(trange[2:T],data_12_pbc_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_12_pbc_std[2:T,:]/sqrt(shots),ylabel="OTOC xx",xlabel="Jt",title="Power Law PBC, h= 12",ylims=[0.,3.2],xlim=[1e-3,1e6])

# ╔═╡ b0ce161e-54ce-4153-abc9-6cf43b49774e
heatmap(1:N,trange[2:T],data_12_pbc_mean[2:T,:],yaxis=:log)

# ╔═╡ 6dd61896-d404-4913-8309-9a9fb1d3121f
md"# NearestNeighbour PBC"

# ╔═╡ d0baa727-8f70-4c03-a899-75387aa0b93c
md"## h = 0"

# ╔═╡ d0c9863d-dbf5-4eff-9cf7-a8ef72ad63d9
begin
	f_0_nn_pbc = "7194615_N13_ED.jld2"
	
	jobids_0_nn_pbc = load(path*f_0_nn_pbc,"jobid")
	params_0_nn_pbc = load(path*f_0_nn_pbc,"params")
	data_0_nn_pbc = 2*ones(T,N,50)-2*load(path*f_0_nn_pbc,"data")

	params_0_nn_pbc
end

# ╔═╡ c459bd9e-d49f-4ab8-8019-4cb10976122c
begin
	data_0_nn_pbc_mean = disorder_mean(data_0_nn_pbc,50)
	data_0_nn_pbc_std = disorder_std(data_0_nn_pbc,50)
end

# ╔═╡ e90ae89c-9331-4610-9fbb-457f9b24d797
plot(trange[2:T],data_0_nn_pbc_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_0_nn_pbc_std[2:T,:]/sqrt(shots),ylabel="OTOC xx",xlabel="Jt",title="Nearest Neighbour PBC, h = 0")

# ╔═╡ cc5b39ed-25c5-4fa5-a53b-9e785e42017d
heatmap(1:N,trange[2:T],data_0_nn_pbc_mean[2:T,:],yaxis=:log)

# ╔═╡ 33ad39fc-1113-4131-ad49-57da9f233d20
md"## h = 4"

# ╔═╡ 4437d9a0-c19c-48a4-b0d9-f52f80737975
begin
	f_4_nn_pbc = "7180027-7194616_N13_ED.jld2"
	combine_files(["7180027_N13_ED.jld2","7194616_N13_ED.jld2"],path,f_4_nn_pbc)
	
	jobids_4_nn_pbc = load(path*f_4_nn_pbc,"jobid")
	params_4_nn_pbc = load(path*f_4_nn_pbc,"params")
	data_4_nn_pbc = 2*ones(T,N,shots)-2*cat(load(path*f_4_nn_pbc,"data")...,dims=3)

	params_4_nn_pbc
end

# ╔═╡ 58867c94-d6e1-4351-8fce-f7790e6ee84e
begin
	data_4_nn_pbc_mean = disorder_mean(data_4_nn_pbc,shots)
	data_4_nn_pbc_std = disorder_std(data_4_nn_pbc,shots)
end

# ╔═╡ 0e7b54d2-4440-4465-a5e4-08552aecc3c9
plot(trange[2:T],data_4_nn_pbc_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_4_nn_pbc_std[2:T,:]/sqrt(shots),ylabel="OTOC xx",xlabel="Jt",title="Nearest Neighbour PBC, h = 4")

# ╔═╡ b502a27e-c36b-47c3-bdd4-f803ba625e57
heatmap(1:N,trange[2:T],data_4_nn_pbc_mean[2:T,:],yaxis=:log)

# ╔═╡ 23945b70-2798-4173-99e9-b434afaae5eb
md"## h = 12"

# ╔═╡ b0a7e802-3f25-4891-b7ff-16b28652bd07
begin
	f_12_nn_pbc = "7180030-7194617_N13_ED.jld2"
	combine_files(["7180030_N13_ED.jld2","7194617_N13_ED.jld2"],path,f_12_nn_pbc)
	
	jobids_12_nn_pbc = load(path*f_12_nn_pbc,"jobid")
	params_12_nn_pbc = load(path*f_12_nn_pbc,"params")
	data_12_nn_pbc = 2*ones(T,N,shots)-2*cat(load(path*f_12_nn_pbc,"data")...,dims=3)

	params_12_nn_pbc
end

# ╔═╡ 3f49b55b-28b9-4306-a836-e5adfe73d950
begin
	data_12_nn_pbc_mean = disorder_mean(data_12_nn_pbc,shots)
	data_12_nn_pbc_std = disorder_std(data_12_nn_pbc,shots)
end

# ╔═╡ d96dce80-c4cc-4ea1-a8b3-916f9609f57d
plot(trange[2:T],data_12_nn_pbc_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_12_nn_pbc_std[2:T,:]/sqrt(shots),ylabel="OTOC xx",xlabel="Jt",title="Nearest Neighbour PBC, h = 12")

# ╔═╡ c237e2de-78b8-434a-86d5-c1df1b289e87
heatmap(1:N,trange[2:T],data_12_nn_pbc_mean[2:T,:],yaxis=:log)

# ╔═╡ b506156f-03f8-43ad-8160-b47542d47ed8
md"# σzσz: N = 13, tmax = 1e10"

# ╔═╡ ad778282-2218-492c-89fe-34ac60e1e341
md"# PowerLaw OBC"

# ╔═╡ 32e658a7-50eb-4207-a11f-ef8ca90ceccc
md"## h = 0"

# ╔═╡ af03157d-5228-4d90-9931-5fc0788ff03b
begin
	f_0zz = "7201192_N13_ED.jld2"
	
	jobids_0zz = load(path*f_0zz,"jobid")
	params_0zz = load(path*f_0zz,"params")
	data_0zz = 2*ones(T,N,50)-2*load(path*f_0zz,"data")

	params_0zz
end

# ╔═╡ d36366a9-4acb-42cb-9943-9e04ad250bb2
begin
	data_0zz_mean = disorder_mean(data_0zz,50)
	data_0zz_std = disorder_std(data_0zz,50)
end

# ╔═╡ b0dc5b3d-e462-4c6d-a7f7-3965a08a8ae9
nearest_neighbourJ_pbc(4)

# ╔═╡ 25cec4e7-c39b-4d55-a200-43ab16ca45a9
plot(trange[2:T],data_0zz_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_0zz_std[2:T,:]/sqrt(50),ylabel="OTOC zz",xlabel="Jt",title="Power Law OBC, h = 0")

# ╔═╡ a0c935e8-a6e4-4fb7-97e6-c1d8200cddea
heatmap(1:N,trange[2:T],data_0zz_mean[2:T,:],yaxis=:log)

# ╔═╡ 1f47b081-7724-4b39-8a71-996f876f2b2c
md"## h = 4"

# ╔═╡ 35c2349b-1e0c-4348-8aca-97401c00e74f
begin
	f_4zz = "7201193-7201194_N13_ED.jld2"
	combine_files(["7201193_N13_ED.jld2","7201194_N13_ED.jld2"],path,f_4zz)
	
	jobids_4zz = load(path*f_4zz,"jobid")
	params_4zz = load(path*f_4zz,"params")
	data_4zz = 2*ones(T,N,shots)-2*cat(load(path*f_4zz,"data")...,dims=3)

	params_4zz
end

# ╔═╡ 763dc4f8-e824-4b88-9dea-417d05c3fbfd
begin
	data_4zz_mean = disorder_mean(data_4zz,shots)
	data_4zz_std = disorder_std(data_4zz,shots)
end

# ╔═╡ 50358394-27e5-4476-a835-b26ba8fc11e5
plot(trange[2:T],data_4zz_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_4zz_std[2:T,:]/sqrt(shots),ylabel="OTOC zz",xlabel="Jt",title="Power Law OBC, h = 4")

# ╔═╡ 6bcf2f5f-d81b-4beb-90a5-73aceae99c08
heatmap(1:N,trange[2:T],data_4zz_mean[2:T,:],yaxis=:log)

# ╔═╡ d715e6ab-858a-452a-997b-cb17b9be3ea8
md"## h = 12"

# ╔═╡ 26e0671c-5fba-4536-8854-a52c3660db73
begin
	f_12zz = "7201195-7201196_N13_ED.jld2"
	combine_files(["7201195_N13_ED.jld2","7201196_N13_ED.jld2"],path,f_12zz)
	
	jobids_12zz = load(path*f_12zz,"jobid")
	params_12zz = load(path*f_12zz,"params")
	data_12zz = 2*ones(T,N,shots)-2*cat(load(path*f_12zz,"data")...,dims=3)

	params_12zz
end

# ╔═╡ 31021ee5-14bc-4b5a-b5bd-12e3d2443478
begin
	data_12zz_mean = disorder_mean(data_12zz,shots)
	data_12zz_std = disorder_std(data_12zz,shots)
end

# ╔═╡ 960bfce8-fc29-4e95-a841-027836125e8d
plot(trange[2:T],data_12zz_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_12zz_std[2:T,:]/sqrt(shots),ylabel="OTOC zz",xlabel="Jt",title="Power Law OBC, h = 12",ylims=[0.,1.2])

# ╔═╡ 9e10cd51-2666-427a-8172-27be5e435f41
heatmap(1:N,trange[2:T],data_12zz_mean[2:T,:],yaxis=:log)

# ╔═╡ c526ce69-027f-478d-bd8b-f69360dc98fc
md"# NearestNeighbour OBC"

# ╔═╡ fa7f51c2-0b07-4841-9425-f3f855ccc081
md"## h = 0"

# ╔═╡ 828c0fea-5559-4382-95b6-761d87c03724
begin
	f_0_nnzz = "7195009_N13_ED.jld2"
	
	jobids_0_nnzz = load(path*f_0_nnzz,"jobid")
	params_0_nnzz = load(path*f_0_nnzz,"params")
	data_0_nnzz = 2*ones(T,N,50)-2*load(path*f_0_nnzz,"data")

	params_0_nnzz
end

# ╔═╡ ecae720a-0623-481c-861a-9be7dc23c152
begin
	data_0_nnzz_mean = disorder_mean(data_0_nnzz,50)
	data_0_nnzz_std = disorder_std(data_0_nnzz,50)
end

# ╔═╡ 8e55c88e-2c64-40c2-aec1-95d3ab423423
plot(trange[2:T],data_0_nnzz_mean[2:110,:],xaxis=:log,legend=:bottomright,ylabel="OTOC zz",xlabel="Jt",title="Nearest Neighbour OBC, h = 0")

# ╔═╡ d9a214fa-a611-4453-a6d1-55196bd722f2
heatmap(1:N,trange[2:T],data_0_nnzz_mean[2:T,:],yaxis=:log)

# ╔═╡ 3aad516e-3bb2-4a9e-ba44-b1f52cafa240
md"## h = 4"

# ╔═╡ 28e679e3-0665-4f54-b346-781fa2487ca5
begin
	f_4_nnzz = "7195010-7195011_N13_ED.jld2"
	combine_files(["7195010_N13_ED.jld2","7195011_N13_ED.jld2"],path,f_4_nnzz)
	
	jobids_4_nnzz = load(path*f_4_nnzz,"jobid")
	params_4_nnzz = load(path*f_4_nnzz,"params")
	data_4_nnzz = 2*ones(T,N,shots)-2*cat(load(path*f_4_nnzz,"data")...,dims=3)

	params_4_nnzz
end

# ╔═╡ 36c36419-b4ab-403c-81d3-0a995487bb7f
begin
	data_4_nnzz_mean = disorder_mean(data_4_nnzz,shots)
	data_4_nnzz_std = disorder_std(data_4_nnzz,shots)
end

# ╔═╡ 9d6cfac8-7221-4bd3-a8a8-2a33e841e63b
plot(trange[2:T],data_4_nnzz_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_4_nnzz_std[2:T,:]/sqrt(shots),ylabel="OTOC zz",xlabel="Jt",title="Nearest Neighbour OBC, h = 4")

# ╔═╡ 2c82a443-12c4-4b5f-846f-26622299363a
heatmap(1:N,trange[2:T],data_4_nnzz_mean[2:T,:],yaxis=:log)

# ╔═╡ 51a6c8af-9953-411d-aaa3-97907c80a50c
md"## h = 12"

# ╔═╡ 45b8a614-4ba9-4063-98a4-367399984dca
begin
	f_12_nnzz = "7195012-7195013_N13_ED.jld2"
	combine_files(["7195012_N13_ED.jld2","7195013_N13_ED.jld2"],path,f_12_nnzz)
	
	jobids_12_nnzz = load(path*f_12_nnzz,"jobid")
	params_12_nnzz= load(path*f_12_nnzz,"params")
	data_12_nnzz = 2*ones(T,N,shots)-2*cat(load(path*f_12_nnzz,"data")...,dims=3)

	params_12_nnzz
end

# ╔═╡ ea7df1cc-e823-4ada-8b6a-378fc2e54f3d
begin
	data_12_nnzz_mean = disorder_mean(data_12_nnzz,shots)
	data_12_nnzz_std = disorder_std(data_12_nnzz,shots)
end

# ╔═╡ 22ff2263-9cde-480c-bdc3-c947022a2acc
plot(trange[2:T],data_12_nnzz_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_12_nnzz_std[2:T,:]/sqrt(shots),ylabel="OTOC zz",xlabel="Jt",title="Nearest Neighbour OBC, h = 12",ylims=[0.,1.2])

# ╔═╡ 36be4821-9712-44d6-a27c-f2e6f6418542
heatmap(1:N,trange[2:T],data_12_nnzz_mean[2:T,:],yaxis=:log)

# ╔═╡ 5c7de876-020e-4ec1-aa03-03bb4cbf4d70
begin
	plot(trange[2:T],data_12_nnzz_mean[2:110,:]-data_12zz_mean[2:110,:],xaxis=:log,legend=:bottomright)#,ribbon=sqrt.(data_12_nnzz_std[2:T,:].^2+data_12zz_std[2:T,:].^2)/sqrt(shots))
end

# ╔═╡ 487a8c45-bed1-45e6-a8d3-39056723597c
md"# PowerLaw PBC"

# ╔═╡ 9370efd2-9b2f-47fd-b3b3-bc1efd0f2dea
md"## h = 0"

# ╔═╡ 2548a2f2-865b-449f-83d4-907089876e43
begin
	f_0_pbczz = "7195091_N13_ED.jld2"
	
	jobids_0_pbczz = load(path*f_0_pbczz,"jobid")
	params_0_pbczz = load(path*f_0_pbczz,"params")
	data_0_pbczz = 2*ones(T,N,50)-2*load(path*f_0_pbczz,"data")

	params_0_pbczz
end

# ╔═╡ 7043d9dc-36be-4bd2-995a-4b42d0e45050
begin
	data_0_pbczz_mean = disorder_mean(data_0_pbczz,50)
	data_0_pbczz_std = disorder_std(data_0_pbczz,50)
end

# ╔═╡ ebcc526b-d747-42f9-9b3a-56dfb37586b2
plot(trange[2:T],data_0_pbczz_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_0_pbczz_std[2:T,:]/sqrt(shots),ylabel="OTOC zz",xlabel="Jt",title="Power Law PBC, h = 0")

# ╔═╡ 34c9f430-95f0-4899-a8e5-53a494582ea0
heatmap(1:N,trange[2:T],data_0_pbczz_mean[2:T,:],yaxis=:log)

# ╔═╡ 3a01e0d2-1bfc-402c-be09-a405983f8620
md"## h = 4"

# ╔═╡ 34dbfa80-13d1-4af0-986a-2456be8f2af4
begin
	f_4_pbczz = "7195111-7195115_N13_ED.jld2"
	combine_files(["7195111_N13_ED.jld2","7195115_N13_ED.jld2"],path,f_4_pbczz)
	
	jobids_4_pbczz = load(path*f_4_pbczz,"jobid")
	params_4_pbczz = load(path*f_4_pbczz,"params")
	data_4_pbczz = 2*ones(T,N,shots)-2*cat(load(path*f_4_pbczz,"data")...,dims=3)

	params_4_pbczz
end

# ╔═╡ 1518c39a-332c-4af3-8a2b-c7593e5235de
begin
	data_4_pbczz_mean = disorder_mean(data_4_pbczz,shots)
	data_4_pbczz_std = disorder_std(data_4_pbczz,shots)
end

# ╔═╡ a6354071-c99c-46aa-a1ba-350f392e9939
plot(trange[2:T],data_4_pbczz_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_4_pbczz_std[2:T,:]/sqrt(shots),ylabel="OTOC zz",xlabel="Jt",title="Power Law PBC, h = 4")

# ╔═╡ 747ec1c5-cd26-474a-9d39-51af16176b4c
heatmap(1:N,trange[2:T],data_4_pbczz_mean[2:T,:],yaxis=:log)

# ╔═╡ 733e5835-9bc8-4149-9ec1-5548bf391d13
md"## h = 12"

# ╔═╡ f7c75d7d-d544-4078-b8e7-4490538c63f6
begin
	f_12_pbczz = "7195116-7195117_N13_ED.jld2"
	combine_files(["7195116_N13_ED.jld2","7195117_N13_ED.jld2"],path,f_12_pbczz)
	
	jobids_12_pbczz = load(path*f_12_pbczz,"jobid")
	params_12_pbczz = load(path*f_12_pbczz,"params")
	data_12_pbczz = 2*ones(T,N,shots)-2*cat(load(path*f_12_pbczz,"data")...,dims=3)

	params_12_pbczz
end

# ╔═╡ 40c6412a-ecc1-450d-ab5d-c58582b4e08a
begin
	data_12_pbczz_mean = disorder_mean(data_12_pbczz,shots)
	data_12_pbczz_std = disorder_std(data_12_pbczz,shots)
end

# ╔═╡ b8d8f62c-56ba-49e4-a0c2-57eac3ee57da
plot(trange[2:T],data_12_pbczz_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_12_pbczz_std[2:T,:]/sqrt(shots),ylabel="OTOC zz",xlabel="Jt",title="Power Law PBC, h = 12",ylims=[0.,1.2])

# ╔═╡ 45def354-c1a9-4e18-a80f-95b375169819
heatmap(1:N,trange[2:T],data_12_pbczz_mean[2:T,:],yaxis=:log)

# ╔═╡ b1b8769e-8c7f-44e7-9992-8cf12c96f090
md"# NearestNeighbour PBC"

# ╔═╡ 88dbb9fa-ee64-497e-9313-43307eb6c688
md"## h = 0"

# ╔═╡ 8abb4021-f5a6-469b-962c-cbcd47038354
begin
	f_0_nn_pbczz = "7195118_N13_ED.jld2"
	
	jobids_0_nn_pbczz = load(path*f_0_nn_pbczz,"jobid")
	params_0_nn_pbczz = load(path*f_0_nn_pbczz,"params")
	data_0_nn_pbczz = 2*ones(T,N,50)-2*load(path*f_0_nn_pbczz,"data")

	params_0_nn_pbczz
end

# ╔═╡ 56cf5180-2849-40ba-8a73-f2bb8a9d0a64
begin
	data_0_nn_pbczz_mean = disorder_mean(data_0_nn_pbczz,50)
	data_0_nn_pbczz_std = disorder_std(data_0_nn_pbczz,50)
end

# ╔═╡ 71d5703f-4eb2-48df-bdb5-931508fb3d42
plot(trange[2:T],data_0_nn_pbczz_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_0_nn_pbczz_std[2:T,:]/sqrt(shots),ylabel="OTOC zz",xlabel="Jt",title="Nearest Neighbour PBC, h = 0")

# ╔═╡ fedae693-83e6-47e2-98b4-d37c7b8cb132
heatmap(1:N,trange[2:T],data_0_nn_pbczz_mean[2:T,:],yaxis=:log)

# ╔═╡ e2c297ed-3a83-44ab-9919-877abe03dcef
md"## h = 4"

# ╔═╡ ab554ab4-3037-4935-9d81-540d54e76804
begin
	f_4_nn_pbczz = "7195119-7195120_N13_ED.jld2"
	combine_files(["7195119_N13_ED.jld2","7195120_N13_ED.jld2"],path,f_4_nn_pbczz)
	
	jobids_4_nn_pbczz = load(path*f_4_nn_pbczz,"jobid")
	params_4_nn_pbczz = load(path*f_4_nn_pbczz,"params")
	data_4_nn_pbczz = 2*ones(T,N,shots)-2*cat(load(path*f_4_nn_pbczz,"data")...,dims=3)

	params_4_nn_pbczz
end

# ╔═╡ ec659f4c-59f1-4e76-8e69-8ff984387817
begin
	data_4_nn_pbczz_mean = disorder_mean(data_4_nn_pbczz,shots)
	data_4_nn_pbczz_std = disorder_std(data_4_nn_pbczz,shots)
end

# ╔═╡ e9d508b3-bb98-4aa6-ad17-bd94dd15f702
plot(trange[2:T],data_4_nn_pbczz_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_4_nn_pbczz_std[2:T,:]/sqrt(shots),ylabel="OTOC zz",xlabel="Jt",title="Nearest Neighbour PBC, h = 4")

# ╔═╡ dbee9aae-a8f1-4dde-a4d7-037933e3d511
heatmap(1:N,trange[2:T],data_4_nn_pbczz_mean[2:T,:],yaxis=:log)

# ╔═╡ 396a7a73-a004-427f-95b7-9090951e18ea
md"## h = 12"

# ╔═╡ 85155eaa-04ca-43b0-a6b4-f027738e2a9f
begin
	f_12_nn_pbczz = "7195121-7195122_N13_ED.jld2"
	combine_files(["7195121_N13_ED.jld2","7195122_N13_ED.jld2"],path,f_12_nn_pbczz)
	
	jobids_12_nn_pbczz = load(path*f_12_nn_pbczz,"jobid")
	params_12_nn_pbczz = load(path*f_12_nn_pbczz,"params")
	data_12_nn_pbczz = 2*ones(T,N,shots)-2*cat(load(path*f_12_nn_pbczz,"data")...,dims=3)

	params_12_nn_pbczz
end

# ╔═╡ 6257ceaf-8956-4de6-b556-a7e3f53f3dec
begin
	data_12_nn_pbczz_mean = disorder_mean(data_12_nn_pbczz,shots)
	data_12_nn_pbczz_std = disorder_std(data_12_nn_pbczz,shots)
end

# ╔═╡ 476f35cb-ed3f-40c7-aec6-efe1ee060607
plot(trange[2:T],data_12_nn_pbczz_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_12_nn_pbczz_std[2:T,:]/sqrt(shots),ylabel="OTOC zz",xlabel="Jt",title="Nearest Neighbour PBC, h = 12")

# ╔═╡ 9ce917fd-3af7-4637-a0cf-1a2a5867a515
heatmap(1:N,trange[2:T],data_12_nn_pbczz_mean[2:T,:],yaxis=:log)

# ╔═╡ Cell order:
# ╠═b630a7e4-3e73-11ed-1baf-c976bb434b06
# ╠═db8bcfe9-3519-43c4-abb8-df5a3a757c73
# ╠═7e135432-4ad6-47c0-a2d6-52da562f66e8
# ╠═166576ac-961d-4a22-9406-327e7e76899d
# ╠═edb5ce15-8d34-4ff3-ab4a-873dbf1a61a3
# ╠═c6561277-446a-49f4-ad4a-a9a60ffe57df
# ╠═3b511bc9-840b-47a7-b4a2-363ca0188e0e
# ╠═ecb93b7f-e7bf-4cfd-b0dd-eab8ecb699fe
# ╠═6b87691e-c010-41c2-ac8d-721c6ce972bb
# ╠═71b9bf79-5135-4ee5-bcdc-7ba1642910af
# ╠═6d7f5c78-8082-491b-a884-5839185641fc
# ╠═6c0d6752-667b-49a7-8f7f-7a039e7104a6
# ╠═4af95c99-f7a3-4d71-b050-a37e8a5f1f03
# ╠═bc844145-f979-4277-a54c-7e902e47ccf0
# ╠═4c422f78-8269-491d-9ba3-2df6bfef198d
# ╠═e22ce722-a048-4df1-b77e-97a113defa53
# ╠═9900588f-d02f-4e46-b136-9903da1c121a
# ╠═88b944a6-ae83-419d-a3b5-988c512bad07
# ╠═fa96cb73-d1cb-401c-addb-d86b8c73dd9a
# ╠═9f5dcb52-8418-417b-b7a3-0f7f7682fcbc
# ╠═a148f020-cb39-4820-b868-503d2dca252c
# ╠═10ca2e1c-d5dc-4c39-bb86-0981d951d283
# ╠═880c16dd-240b-4ad5-92f1-54b923ef75da
# ╠═8bbe06ac-a2da-4ceb-bc5c-3c442a78c54c
# ╠═be7d456d-b84d-47fd-a750-949227cbea71
# ╠═e8d0bd88-11aa-419a-98bc-05c3cfc3cb13
# ╠═c4f385eb-6e67-402d-9f98-c8ac8621c9a2
# ╠═4302e49d-d06a-4e22-86a1-75669096779f
# ╠═9380deb5-848b-40f9-a58d-0a8d157c90b1
# ╠═6602f0c5-1445-417a-b2ca-490622c03169
# ╠═bcc17df1-370a-4151-9433-688f3195bc8e
# ╠═6f8b67b3-796f-4246-a250-0a023f517ace
# ╠═a939797d-8e74-4d0f-bfff-deec14b7ac51
# ╠═33282538-dadb-4487-986e-ea89f33ebe88
# ╠═ac5908dc-61e9-479b-a33c-50e70abd829f
# ╠═72391c41-6e9a-4c4c-a766-57b27d10445e
# ╠═8679c844-ceb7-44c2-9e88-e13456c77933
# ╠═0d72b7be-8a28-400d-9c55-a9ada9d62322
# ╠═e7af56ef-24ed-4e49-917e-b05fa0473c0c
# ╠═a573b238-1a93-46d1-b7f5-3deccc02bdc5
# ╠═b72312cb-98a8-4a93-bd9f-8a5b5e8c6750
# ╠═88300d49-ecb9-4f08-9262-f109b8297b03
# ╠═44181015-34a4-47ef-aefb-6f8509fdfca1
# ╠═50f76e36-ec8b-4989-bfe2-d109f811ae26
# ╠═54d6b8bc-6b9b-4243-9aea-56eafbf9f701
# ╠═b969969c-967e-4ff6-b781-62087e167f29
# ╠═d21541e3-ba62-4dd5-844e-1751cf134d29
# ╠═f55a5ebe-b8cf-4772-8aec-ee201db139ef
# ╠═7864764e-8f34-43b9-bedd-970847d38688
# ╠═c9d599c4-6cd8-4733-9603-8bfffd3ac836
# ╠═9a75ae46-c32d-4d0d-a7d6-e1df83f17129
# ╠═9dcfe111-e17f-4b7c-807b-a5692d88083a
# ╠═0a04bf3f-9783-4d1d-99e3-a00e7d844d07
# ╠═5f376697-c573-4871-b8d1-758f9b46b92b
# ╠═386d9b58-bf17-49d9-a2b8-443aad198c6c
# ╠═675bca7e-d84a-4774-89f8-9f67ea3e4d2c
# ╠═cbbe5b27-2a57-4296-a4ef-b87c0502af4d
# ╠═da1670a7-070c-4533-bef5-60b07a3f5edf
# ╠═d7e301c3-4e11-444e-80ac-8691b9d4c556
# ╠═727e60b7-ec48-4a13-b4c7-48b68f27f9e6
# ╠═96da9674-935c-444a-81ae-a98fa54a95d4
# ╠═b2af0bf9-f305-48c6-8295-c560bb9d25ec
# ╠═9378d2bd-2e31-4a7d-a81c-574e70c9a171
# ╠═b0ce161e-54ce-4153-abc9-6cf43b49774e
# ╠═6dd61896-d404-4913-8309-9a9fb1d3121f
# ╠═d0baa727-8f70-4c03-a899-75387aa0b93c
# ╠═d0c9863d-dbf5-4eff-9cf7-a8ef72ad63d9
# ╠═c459bd9e-d49f-4ab8-8019-4cb10976122c
# ╠═e90ae89c-9331-4610-9fbb-457f9b24d797
# ╠═cc5b39ed-25c5-4fa5-a53b-9e785e42017d
# ╠═33ad39fc-1113-4131-ad49-57da9f233d20
# ╠═4437d9a0-c19c-48a4-b0d9-f52f80737975
# ╠═58867c94-d6e1-4351-8fce-f7790e6ee84e
# ╠═0e7b54d2-4440-4465-a5e4-08552aecc3c9
# ╠═b502a27e-c36b-47c3-bdd4-f803ba625e57
# ╠═23945b70-2798-4173-99e9-b434afaae5eb
# ╠═b0a7e802-3f25-4891-b7ff-16b28652bd07
# ╠═3f49b55b-28b9-4306-a836-e5adfe73d950
# ╠═d96dce80-c4cc-4ea1-a8b3-916f9609f57d
# ╠═c237e2de-78b8-434a-86d5-c1df1b289e87
# ╠═b506156f-03f8-43ad-8160-b47542d47ed8
# ╠═ad778282-2218-492c-89fe-34ac60e1e341
# ╠═32e658a7-50eb-4207-a11f-ef8ca90ceccc
# ╠═af03157d-5228-4d90-9931-5fc0788ff03b
# ╠═d36366a9-4acb-42cb-9943-9e04ad250bb2
# ╠═b0dc5b3d-e462-4c6d-a7f7-3965a08a8ae9
# ╠═25cec4e7-c39b-4d55-a200-43ab16ca45a9
# ╠═a0c935e8-a6e4-4fb7-97e6-c1d8200cddea
# ╠═1f47b081-7724-4b39-8a71-996f876f2b2c
# ╠═35c2349b-1e0c-4348-8aca-97401c00e74f
# ╠═763dc4f8-e824-4b88-9dea-417d05c3fbfd
# ╠═50358394-27e5-4476-a835-b26ba8fc11e5
# ╠═6bcf2f5f-d81b-4beb-90a5-73aceae99c08
# ╠═d715e6ab-858a-452a-997b-cb17b9be3ea8
# ╠═26e0671c-5fba-4536-8854-a52c3660db73
# ╠═31021ee5-14bc-4b5a-b5bd-12e3d2443478
# ╠═960bfce8-fc29-4e95-a841-027836125e8d
# ╠═9e10cd51-2666-427a-8172-27be5e435f41
# ╠═c526ce69-027f-478d-bd8b-f69360dc98fc
# ╠═fa7f51c2-0b07-4841-9425-f3f855ccc081
# ╠═828c0fea-5559-4382-95b6-761d87c03724
# ╠═ecae720a-0623-481c-861a-9be7dc23c152
# ╠═8e55c88e-2c64-40c2-aec1-95d3ab423423
# ╠═d9a214fa-a611-4453-a6d1-55196bd722f2
# ╠═3aad516e-3bb2-4a9e-ba44-b1f52cafa240
# ╠═28e679e3-0665-4f54-b346-781fa2487ca5
# ╠═36c36419-b4ab-403c-81d3-0a995487bb7f
# ╠═9d6cfac8-7221-4bd3-a8a8-2a33e841e63b
# ╠═2c82a443-12c4-4b5f-846f-26622299363a
# ╠═51a6c8af-9953-411d-aaa3-97907c80a50c
# ╠═45b8a614-4ba9-4063-98a4-367399984dca
# ╠═ea7df1cc-e823-4ada-8b6a-378fc2e54f3d
# ╠═22ff2263-9cde-480c-bdc3-c947022a2acc
# ╠═36be4821-9712-44d6-a27c-f2e6f6418542
# ╠═5c7de876-020e-4ec1-aa03-03bb4cbf4d70
# ╠═487a8c45-bed1-45e6-a8d3-39056723597c
# ╠═9370efd2-9b2f-47fd-b3b3-bc1efd0f2dea
# ╠═2548a2f2-865b-449f-83d4-907089876e43
# ╠═7043d9dc-36be-4bd2-995a-4b42d0e45050
# ╠═ebcc526b-d747-42f9-9b3a-56dfb37586b2
# ╠═34c9f430-95f0-4899-a8e5-53a494582ea0
# ╠═3a01e0d2-1bfc-402c-be09-a405983f8620
# ╠═34dbfa80-13d1-4af0-986a-2456be8f2af4
# ╠═1518c39a-332c-4af3-8a2b-c7593e5235de
# ╠═a6354071-c99c-46aa-a1ba-350f392e9939
# ╠═747ec1c5-cd26-474a-9d39-51af16176b4c
# ╠═733e5835-9bc8-4149-9ec1-5548bf391d13
# ╠═f7c75d7d-d544-4078-b8e7-4490538c63f6
# ╠═40c6412a-ecc1-450d-ab5d-c58582b4e08a
# ╠═b8d8f62c-56ba-49e4-a0c2-57eac3ee57da
# ╠═45def354-c1a9-4e18-a80f-95b375169819
# ╠═b1b8769e-8c7f-44e7-9992-8cf12c96f090
# ╠═88dbb9fa-ee64-497e-9313-43307eb6c688
# ╠═8abb4021-f5a6-469b-962c-cbcd47038354
# ╠═56cf5180-2849-40ba-8a73-f2bb8a9d0a64
# ╠═71d5703f-4eb2-48df-bdb5-931508fb3d42
# ╠═fedae693-83e6-47e2-98b4-d37c7b8cb132
# ╠═e2c297ed-3a83-44ab-9919-877abe03dcef
# ╠═ab554ab4-3037-4935-9d81-540d54e76804
# ╠═ec659f4c-59f1-4e76-8e69-8ff984387817
# ╠═e9d508b3-bb98-4aa6-ad17-bd94dd15f702
# ╠═dbee9aae-a8f1-4dde-a4d7-037933e3d511
# ╠═396a7a73-a004-427f-95b7-9090951e18ea
# ╠═85155eaa-04ca-43b0-a6b4-f027738e2a9f
# ╠═6257ceaf-8956-4de6-b556-a7e3f53f3dec
# ╠═476f35cb-ed3f-40c7-aec6-efe1ee060607
# ╠═9ce917fd-3af7-4637-a0cf-1a2a5867a515
