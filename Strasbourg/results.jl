### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ c696705d-73cc-4d72-a51c-4c51446b9d64
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 92fb6306-4024-11ed-1377-5bd63b6ff1ca
using LinearAlgebra,Plots,JLD2,Statistics,PlutoUI

# ╔═╡ d09eb176-872e-4a05-8625-e1b7befb9e81
TableOfContents()

# ╔═╡ 57ed6799-fb4c-452b-9e99-71640d7692cb
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

# ╔═╡ 3a797afb-78fd-47d0-9c1c-fb3ef92a098c
begin
	N = 12
	path = pwd() *"/data/"
	shots = 100
	trange = logrange(-2,10,1e10)
	T = length(trange)
end

# ╔═╡ 43a24fc3-ae17-4a08-a22f-15011da7b410
md"# σxσx"

# ╔═╡ 4f25a173-0fe3-4e70-9e9f-df2d6a7259b4
md"## XXZ"

# ╔═╡ 24e973f7-d493-41df-ac01-38863668c0de
begin
	f_xxz = "7205359_N12_ED.jld2"
	f_xxzneel = "7178194_N12_ED.jld2"
	f_xxzneelx = "7178196_N12_ED.jld2"
	
	jobids_xxz = load(path*f_xxz,"jobid")
	params_xxz = load(path*f_xxz,"params")
	data_xxz = 2*ones(T,N,shots)-2*load(path*f_xxz,"data")

	jobids_xxzneel = load(path*f_xxzneel,"jobid")
	params_xxzneel = load(path*f_xxzneel,"params")
	data_xxzneel = 2*ones(T,N,shots)-2*load(path*f_xxzneel,"data")

	jobids_xxzneelx = load(path*f_xxzneelx,"jobid")
	params_xxzneelx = load(path*f_xxzneelx,"params")
	data_xxzneelx = 2*ones(T,N,shots)-2*load(path*f_xxzneelx,"data")

	size(data_xxz)
end

# ╔═╡ 7eb921df-184f-40b1-839c-5ccca39b0fa2
begin
	data_xxz_mean = disorder_mean(data_xxz,shots)
	data_xxz_std = disorder_std(data_xxz,shots)

	data_xxzneel_mean = disorder_mean(data_xxzneel,shots)
	data_xxzneel_std = disorder_std(data_xxzneel,shots)

	data_xxzneelx_mean = disorder_mean(data_xxzneelx,shots)
	data_xxzneelx_std = disorder_std(data_xxzneelx,shots)
end

# ╔═╡ c63bfcbd-9369-416f-8720-83efd2873cf9
plot(trange[2:T],data_xxz_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_xxz_std[2:T,:]/sqrt(shots),ylabel="OTOC xx",xlabel="Jt",title="XXZ Nearest Neighbour OBC, h = 12")

# ╔═╡ 4cdf7961-1938-4db1-a1ab-35d9e96ffc0a
plot(trange[2:T],data_xxzneel_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_xxzneel_std[2:T,:]/sqrt(shots),ylabel="OTOC xx",xlabel="Jt",title="XXZ Nearest Neighbour OBC, h = 12, Néel")

# ╔═╡ dd971d85-b6e1-4772-a842-da7010d90d62
plot(trange[2:T],data_xxzneelx_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_xxzneelx_std[2:T,:]/sqrt(shots),ylabel="OTOC xx",xlabel="Jt",title="XXZ Nearest Neighbour OBC, h = 12, Néel x")

# ╔═╡ 8520d920-c4c6-409e-9ea3-92279e6ef203
heatmap(1:N,trange[2:T],data_xxz_mean[2:T,:],yaxis=:log)

# ╔═╡ 20ceb700-dc96-4102-8e47-beea3882974f
heatmap(1:N,trange[2:T],data_xxzneel_mean[2:T,:],yaxis=:log)

# ╔═╡ daed58b3-ce9c-4ecd-93fc-57935fcad5f4
heatmap(1:N,trange[2:T],data_xxzneelx_mean[2:T,:],yaxis=:log)

# ╔═╡ 3faaa2d3-d06e-46a7-8fc4-3380d6835639
md"## XYZ"

# ╔═╡ 8ef450cb-5321-4765-ad62-4c2e01b51c8a
begin
	f_xyz = "7205360_N12_ED.jld2"
	f_xyzneel = "7178197_N12_ED.jld2"
	
	jobids_xyz = load(path*f_xyz,"jobid")
	params_xyz = load(path*f_xyz,"params")
	data_xyz = 2*ones(T,N,shots)-2*load(path*f_xyz,"data")

	jobids_xyzneel = load(path*f_xyzneel,"jobid")
	params_xyzneel = load(path*f_xyzneel,"params")
	data_xyzneel = 2*ones(T,N,shots)-2*load(path*f_xyzneel,"data")

	size(data_xyz)
end

# ╔═╡ dec32486-ccba-436c-9b87-aa62bbb82dc5
begin
	data_xyz_mean = disorder_mean(data_xyz,shots)
	data_xyz_std = disorder_std(data_xyz,shots)

	data_xyzneel_mean = disorder_mean(data_xyzneel,shots)
	data_xyzneel_std = disorder_std(data_xyzneel,shots)
end

# ╔═╡ 620b1a27-cc10-4480-90da-23a85044f56c
plot(trange[2:T],data_xyz_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_xyz_std[2:T,:]/sqrt(shots),ylabel="OTOC xx",xlabel="Jt",title="XYZ Nearest Neighbour OBC, h = 12")

# ╔═╡ 5ea785f8-00f0-46f5-8e55-6d14847c6877
plot(trange[2:T],data_xyzneel_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_xyzneel_std[2:T,:]/sqrt(shots),ylabel="OTOC xx",xlabel="Jt",title="XYZ Nearest Neighbour OBC, h = 12, Néel")

# ╔═╡ 0945e212-9bdd-48c8-89ca-cd1b32cbd4c9
heatmap(1:N,trange[2:T],data_xyz_mean[2:T,:],yaxis=:log)

# ╔═╡ c8132ff6-e867-4c5c-8f07-ab96ce40664a
heatmap(1:N,trange[2:T],data_xyzneel_mean[2:T,:],yaxis=:log)

# ╔═╡ 1fc9afda-c500-4a74-bb8e-3800c5f2e486
md"## TFIM"

# ╔═╡ 5f147280-7f52-480d-a55a-4e1b828264fb
begin
	f_tfim = "7205361_N12_ED.jld2"
	f_tfimneel = "7178521_N12_ED.jld2"
	
	jobids_tfim = load(path*f_tfim,"jobid")
	params_tfim = load(path*f_tfim,"params")
	data_tfim = 2*ones(T,N,shots)-2*load(path*f_tfim,"data")

	jobids_tfimneel = load(path*f_tfimneel,"jobid")
	params_tfimneel = load(path*f_tfimneel,"params")
	data_tfimneel = 2*ones(T,N,shots)-2*load(path*f_tfimneel,"data")

	size(data_tfim)
end

# ╔═╡ 02d82c05-1cc2-4eeb-b04a-9792f5c695bc
begin
	data_tfim_mean = disorder_mean(data_tfim,shots)
	data_tfim_std = disorder_std(data_tfim,shots)

	data_tfimneel_mean = disorder_mean(data_tfim,shots)
	data_tfimneel_std = disorder_std(data_tfim,shots)
end

# ╔═╡ 24b68887-e80d-4ee2-9670-91660b941522
plot(trange[2:T],data_tfim_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_tfim_std[2:T,:]/sqrt(shots),ylabel="OTOC xx",xlabel="Jt",title="TFIM Nearest Neighbour OBC, h = 12")

# ╔═╡ 32ea0426-aa25-45f6-920a-515bb97376e3
plot(trange[2:T],data_tfimneel_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_tfimneel_std[2:T,:]/sqrt(shots),ylabel="OTOC xx",xlabel="Jt",title="TFIM Nearest Neighbour OBC, h = 12, Néel")

# ╔═╡ 751f7ed5-45de-4174-8476-c427a96192ca
heatmap(1:N,trange[2:T],data_tfim_mean[2:T,:],yaxis=:log)

# ╔═╡ 557fc84c-fb9b-4f93-8a5f-af7595eab5e1
heatmap(1:N,trange[2:T],data_tfimneel_mean[2:T,:],yaxis=:log)

# ╔═╡ 63884925-bf8f-4e1d-b414-1539a84f345d
md"# σzσz"

# ╔═╡ b78a012d-8ee3-4e62-a5af-2582f0a05d20
md"## XXZ"

# ╔═╡ 36b03fd2-4d92-4e0e-a6a4-ccdd5aa6c4bc
begin
	f_xxz2 = "7216358_N12_ED.jld2"
	f_xxzneel2 = "7216361_N12_ED.jld2"
	f_xxzneelx2 = "7216364_N12_ED.jld2"
	
	jobids_xxz2 = load(path*f_xxz2,"jobid")
	params_xxz2 = load(path*f_xxz2,"params")
	data_xxz2 = 2*ones(T,N,shots)-2*load(path*f_xxz2,"data")

	jobids_xxzneel2 = load(path*f_xxzneel2,"jobid")
	params_xxzneel2 = load(path*f_xxzneel2,"params")
	data_xxzneel2 = 2*ones(T,N,shots)-2*load(path*f_xxzneel2,"data")

	jobids_xxzneelx2 = load(path*f_xxzneelx2,"jobid")
	params_xxzneelx2 = load(path*f_xxzneelx2,"params")
	data_xxzneelx2 = 2*ones(T,N,shots)-2*load(path*f_xxzneelx2,"data")

	size(data_xxz2)
end

# ╔═╡ ed16066e-f998-4cc7-9987-b6862ea7140a
begin
	data_xxz2_mean = disorder_mean(data_xxz2,shots)
	data_xxz2_std = disorder_std(data_xxz2,shots)

	data_xxzneel2_mean = disorder_mean(data_xxzneel2,shots)
	data_xxzneel2_std = disorder_std(data_xxzneel2,shots)

	data_xxzneelx2_mean = disorder_mean(data_xxzneelx2,shots)
	data_xxzneelx2_std = disorder_std(data_xxzneelx2,shots)
end

# ╔═╡ 0326d95e-9d1a-4059-a558-77440b9a6462
plot(trange[2:T],data_xxz2_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_xxz2_std[2:T,:]/sqrt(shots),ylabel="OTOC zz",xlabel="Jt",title="XXZ Nearest Neighbour OBC, h = 12")

# ╔═╡ ffdf7217-6ca7-4f1d-bffc-5ea53fe05dbb
plot(trange[2:T],data_xxzneel2_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_xxzneel2_std[2:T,:]/sqrt(shots),ylabel="OTOC zz",xlabel="Jt",title="XXZ Nearest Neighbour OBC, h = 12, Néel")

# ╔═╡ 7db7dfdc-c1ed-4fb4-8845-79d5766280c6
plot(trange[2:T],data_xxzneelx2_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_xxzneelx2_std[2:T,:]/sqrt(shots),ylabel="OTOC zz",xlabel="Jt",title="XXZ Nearest Neighbour OBC, h = 12, Néel x")

# ╔═╡ 235afcc0-6448-42e5-bc1a-e4efd3a641d8
heatmap(1:N,trange[2:T],data_xxz2_mean[2:T,:],yaxis=:log)

# ╔═╡ dc955f29-3c51-4227-8315-1a611a1eadf5
heatmap(1:N,trange[2:T],data_xxzneel2_mean[2:T,:],yaxis=:log)

# ╔═╡ ea4bab68-56e5-4e5c-8730-f88c601d648e
heatmap(1:N,trange[2:T],data_xxzneelx2_mean[2:T,:],yaxis=:log)

# ╔═╡ 5ed520fa-b147-4aa9-84af-4470f28022a9
md"## XYZ"

# ╔═╡ 232fabb2-da4a-4d15-9f71-8b54c193a915
begin
	f_xyz2 = "7216359_N12_ED.jld2"
	f_xyzneel2 = "7216362_N12_ED.jld2"
	
	jobids_xyz2 = load(path*f_xyz2,"jobid")
	params_xyz2 = load(path*f_xyz2,"params")
	data_xyz2 = 2*ones(T,N,shots)-2*load(path*f_xyz2,"data")

	jobids_xyzneel2 = load(path*f_xyzneel2,"jobid")
	params_xyzneel2 = load(path*f_xyzneel2,"params")
	data_xyzneel2 = 2*ones(T,N,shots)-2*load(path*f_xyzneel2,"data")

	size(data_xyz2)
end

# ╔═╡ 2f7081d9-92cd-454c-abe3-eb795a713488
begin
	data_xyz2_mean = disorder_mean(data_xyz2,shots)
	data_xyz2_std = disorder_std(data_xyz2,shots)

	data_xyzneel2_mean = disorder_mean(data_xyzneel2,shots)
	data_xyzneel2_std = disorder_std(data_xyzneel2,shots)
end

# ╔═╡ 68685054-36a2-4e42-bcc0-6891ccc1a194
plot(trange[2:T],data_xyz2_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_xyz2_std[2:T,:]/sqrt(shots),ylabel="OTOC zz",xlabel="Jt",title="XYZ Nearest Neighbour OBC, h = 12")

# ╔═╡ 2cbe478f-1515-4409-aab9-e258039ea6e5
plot(trange[2:T],data_xyzneel2_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_xyzneel2_std[2:T,:]/sqrt(shots),ylabel="OTOC zz",xlabel="Jt",title="XYZ Nearest Neighbour OBC, h = 12, Néel")

# ╔═╡ 5cd19f95-1be0-4a58-a41b-746b09fb65f6
heatmap(1:N,trange[2:T],data_xyz2_mean[2:T,:],yaxis=:log)

# ╔═╡ b5476531-dc5e-4eb4-9b4b-176aa454e3f3
heatmap(1:N,trange[2:T],data_xyzneel2_mean[2:T,:],yaxis=:log)

# ╔═╡ baf78f11-3e9f-4285-a644-1ef7b8793338
md"## TFIM"

# ╔═╡ 281bdec5-0f2c-4533-8c5e-bd338592a6e8
begin
	f_tfim2 = "7216360_N12_ED.jld2"
	f_tfimneel2 = "7216395_N12_ED.jld2"
	
	jobids_tfim2 = load(path*f_tfim2,"jobid")
	params_tfim2 = load(path*f_tfim2,"params")
	data_tfim2 = 2*ones(T,N,shots)-2*load(path*f_tfim2,"data")

	jobids_tfimneel2 = load(path*f_tfimneel2,"jobid")
	params_tfimneel2 = load(path*f_tfimneel2,"params")
	data_tfimneel2 = 2*ones(T,N,shots)-2*load(path*f_tfimneel2,"data")

	size(data_tfim2)
end

# ╔═╡ 09fd8d17-8d0c-4e76-bfa1-a60695e405f8
begin
	data_tfim2_mean = disorder_mean(data_tfim2,shots)
	data_tfim2_std = disorder_std(data_tfim2,shots)

	data_tfimneel2_mean = disorder_mean(data_tfim2,shots)
	data_tfimneel2_std = disorder_std(data_tfim2,shots)
end

# ╔═╡ 05504117-0cc7-4769-b2bb-ea398987a659
plot(trange[2:T],data_tfim2_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_tfim2_std[2:T,:]/sqrt(shots),ylabel="OTOC zz",xlabel="Jt",title="TFIM Nearest Neighbour OBC, h = 12")

# ╔═╡ 3877d57e-e394-4fdb-98ab-8bd5fd86eaa1
plot(trange[2:T],data_tfimneel2_mean[2:110,:],xaxis=:log,legend=:bottomright,ribbon=data_tfimneel2_std[2:T,:]/sqrt(shots),ylabel="OTOC xx",xlabel="Jt",title="TFIM Nearest Neighbour OBC, h = 12, Néel")

# ╔═╡ 812e9ef6-5450-4228-90be-a5047861d5fd
heatmap(1:N,trange[2:T],data_tfim2_mean[2:T,:],yaxis=:log)

# ╔═╡ 0e026f0f-fc71-4744-8525-7603cb3642f9
heatmap(1:N,trange[2:T],data_tfimneel2_mean[2:T,:],yaxis=:log)

# ╔═╡ 5a8a968f-c410-4017-989a-b2a2c2c8323f


# ╔═╡ 3b7bc143-052f-43a5-b10d-ec622de83d42


# ╔═╡ 314a2063-7580-4453-abfd-32c52799fa34


# ╔═╡ Cell order:
# ╠═92fb6306-4024-11ed-1377-5bd63b6ff1ca
# ╠═c696705d-73cc-4d72-a51c-4c51446b9d64
# ╠═d09eb176-872e-4a05-8625-e1b7befb9e81
# ╠═57ed6799-fb4c-452b-9e99-71640d7692cb
# ╠═3a797afb-78fd-47d0-9c1c-fb3ef92a098c
# ╠═43a24fc3-ae17-4a08-a22f-15011da7b410
# ╠═4f25a173-0fe3-4e70-9e9f-df2d6a7259b4
# ╠═24e973f7-d493-41df-ac01-38863668c0de
# ╠═7eb921df-184f-40b1-839c-5ccca39b0fa2
# ╠═c63bfcbd-9369-416f-8720-83efd2873cf9
# ╠═4cdf7961-1938-4db1-a1ab-35d9e96ffc0a
# ╠═dd971d85-b6e1-4772-a842-da7010d90d62
# ╠═8520d920-c4c6-409e-9ea3-92279e6ef203
# ╠═20ceb700-dc96-4102-8e47-beea3882974f
# ╠═daed58b3-ce9c-4ecd-93fc-57935fcad5f4
# ╠═3faaa2d3-d06e-46a7-8fc4-3380d6835639
# ╠═8ef450cb-5321-4765-ad62-4c2e01b51c8a
# ╠═dec32486-ccba-436c-9b87-aa62bbb82dc5
# ╠═620b1a27-cc10-4480-90da-23a85044f56c
# ╠═5ea785f8-00f0-46f5-8e55-6d14847c6877
# ╠═0945e212-9bdd-48c8-89ca-cd1b32cbd4c9
# ╠═c8132ff6-e867-4c5c-8f07-ab96ce40664a
# ╠═1fc9afda-c500-4a74-bb8e-3800c5f2e486
# ╠═5f147280-7f52-480d-a55a-4e1b828264fb
# ╠═02d82c05-1cc2-4eeb-b04a-9792f5c695bc
# ╠═24b68887-e80d-4ee2-9670-91660b941522
# ╠═32ea0426-aa25-45f6-920a-515bb97376e3
# ╠═751f7ed5-45de-4174-8476-c427a96192ca
# ╠═557fc84c-fb9b-4f93-8a5f-af7595eab5e1
# ╠═63884925-bf8f-4e1d-b414-1539a84f345d
# ╠═b78a012d-8ee3-4e62-a5af-2582f0a05d20
# ╠═36b03fd2-4d92-4e0e-a6a4-ccdd5aa6c4bc
# ╠═ed16066e-f998-4cc7-9987-b6862ea7140a
# ╠═0326d95e-9d1a-4059-a558-77440b9a6462
# ╠═ffdf7217-6ca7-4f1d-bffc-5ea53fe05dbb
# ╠═7db7dfdc-c1ed-4fb4-8845-79d5766280c6
# ╠═235afcc0-6448-42e5-bc1a-e4efd3a641d8
# ╠═dc955f29-3c51-4227-8315-1a611a1eadf5
# ╠═ea4bab68-56e5-4e5c-8730-f88c601d648e
# ╠═5ed520fa-b147-4aa9-84af-4470f28022a9
# ╠═232fabb2-da4a-4d15-9f71-8b54c193a915
# ╠═2f7081d9-92cd-454c-abe3-eb795a713488
# ╠═68685054-36a2-4e42-bcc0-6891ccc1a194
# ╠═2cbe478f-1515-4409-aab9-e258039ea6e5
# ╠═5cd19f95-1be0-4a58-a41b-746b09fb65f6
# ╠═b5476531-dc5e-4eb4-9b4b-176aa454e3f3
# ╠═baf78f11-3e9f-4285-a644-1ef7b8793338
# ╠═281bdec5-0f2c-4533-8c5e-bd338592a6e8
# ╠═09fd8d17-8d0c-4e76-bfa1-a60695e405f8
# ╠═05504117-0cc7-4769-b2bb-ea398987a659
# ╠═3877d57e-e394-4fdb-98ab-8bd5fd86eaa1
# ╠═812e9ef6-5450-4228-90be-a5047861d5fd
# ╠═0e026f0f-fc71-4744-8525-7603cb3642f9
# ╠═5a8a968f-c410-4017-989a-b2a2c2c8323f
# ╠═3b7bc143-052f-43a5-b10d-ec622de83d42
# ╠═314a2063-7580-4453-abfd-32c52799fa34
