### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ c39afa38-c3d8-49cd-b730-3919e7b7a069
using LinearAlgebra,Plots,JLD2,Statistics,PlutoUI

# ╔═╡ 18a0434c-bacf-411f-a840-fb69a9a91c7c
md"## Random State"

# ╔═╡ 6aec44d6-1e21-4f5d-85f8-a25cd942e0c4
pwd()*"/RS"

# ╔═╡ 6612a15a-c3ca-49da-823f-808e8fc6f8ac
begin
	path_RS = pwd()*"/RS/"
	F_RS = 13
	jobids_RS = Vector{String}(undef,F_RS)
	params_RS = Vector{Any}(undef,F_RS)
	data_RS = Vector{Array{Float64,4}}(undef,F_RS)
	N_RS = [8,9,10,11,12,13,14,15,16,17,18,19,20]
	files_RS = ["5213301_N8.jld2","5213366_N9.jld2","5179026_N10.jld2","5179027_N11.jld2","5179028_N12.jld2","5179029_N13.jld2","5173899_N14.jld2","5173900_N15.jld2","5173901_N16.jld2","5173902_N17.jld2","5173903_N18.jld2","5173904_N19.jld2","5173905_N20.jld2"]
	for (i,f) in enumerate(files_RS)
		jobids_RS[i] = load(path_RS*f,"jobid")
		params_RS[i] = load(path_RS*f,"params")
		data_RS[i] = 2*ones(51,N_RS[i],1,10)-2*load(path_RS*f,"data")
	end
end

# ╔═╡ 2e9d9a87-b8ee-48f4-bfcd-33514a9c0250
md"## Larger Sample Random State"

# ╔═╡ e662f523-7129-4ea7-8564-ab2473449b7c
begin
	Test = ["1","2","3"]
	Test[1:length(Test)-1]
end

# ╔═╡ 6f6f4f0f-d83d-48c8-afed-9b2ba5c02523
[x for x in Test]

# ╔═╡ a52bd686-43c8-49af-993d-174a5d10119b
begin
	path_L = pwd()*"/100S/"
	F_L = 7
	jobids_L = Vector{String}(undef,F_L)
	params_L = Vector{Any}(undef,F_L)
	data_L = Vector{Array{Float64,4}}(undef,F_L)
	N_L = [8,9,10,11,12,13,14]
	files_L = ["5203834_N8.jld2","5203835_N9.jld2","5203836_N10.jld2","5203837_N11.jld2","5203838_N12.jld2","5203839_N13.jld2","5203840_N14.jld2"]
	for (i,f) in enumerate(files_L)
		jobids_L[i] = load(path_L*f,"jobid")
		params_L[i] = load(path_L*f,"params")
		data_L[i] = 2*ones(51,N_L[i],1,100)-2*load(path_L*f,"data")
	end
end

# ╔═╡ 79d282be-08d3-4bf3-b700-eee967e1815a
md"## Random Product State"

# ╔═╡ f636cf34-b778-4259-abd1-bede4cf91a0f


# ╔═╡ 64ac6e06-1c34-4430-9bc4-d0bb9097834f
function combine_files(files,path,new_file)
	jobids = Vector{String}(undef,length(files))
	params = Vector{Any}(undef,length(files))
	data = Vector{Array{Float64,4}}(undef,length(files))
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

# ╔═╡ 397f89e9-2051-4109-90ac-47e31a993637
begin
	path_RPS = pwd()*"/RPS/"
	F_RPS = 13
	jobids_RPS = Vector{Any}(undef,F_RPS)
	params_RPS = Vector{Any}(undef,F_RPS)
	data_RPS = Vector{Array{Float64,4}}(undef,F_RPS)
	N_RPS = [8,9,10,11,12,13,14,15,16,17,18,19,20]
	combine_files(["5213440_N17.jld2","5213441_N17.jld2","5213442_N17.jld2"],path_RPS,"5213440-5213442.jld2")
	combine_files(["5213444_N18.jld2","5213445_N18.jld2","5213446_N18.jld2","5213447_N18.jld2","5213448_N18.jld2"],path_RPS,"5213444-5213448.jld2")
	combine_files(["5213531_N19.jld2","5213535_N19.jld2","5213544_N19.jld2","5213552_N19.jld2","5213553_N19.jld2"],path_RPS,"5213531-5213553.jld2")
	combine_files(["5213555_N20.jld2","5213556_N20.jld2","5213557_N20.jld2","5213558_N20.jld2","5213559_N20.jld2","5213560_N20.jld2","5213561_N20.jld2","5213562_N20.jld2","5213563_N20.jld2","5213564_N20.jld2"],path_RPS,"5213555-5213564.jld2")
	files_RPS = ["5203820_N8.jld2","5203822_N9.jld2","5203823_N10.jld2","5203824_N11.jld2","5203825_N12.jld2","5203826_N13.jld2","5203827_N14.jld2","5213368_N15.jld2","5213371_N16.jld2","5213440-5213442.jld2","5213444-5213448.jld2","5213531-5213553.jld2","5213555-5213564.jld2"]
	for (i,f) in enumerate(files_RPS)
		jobids_RPS[i] = load(path_RPS*f,"jobid")
		params_RPS[i] = load(path_RPS*f,"params")
		if N_RPS[i] <= 16
			data_RPS[i] = 2*ones(51,N_RPS[i],1,100)-2*load(path_RPS*f,"data")
		elseif N_RPS[i] >= 17
			data_RPS[i] = 2*ones(51,N_RPS[i],1,100)-2*cat(load(path_RPS*f,"data")...,dims=4)[:,:,:,1:100]
		else
			data_RPS[i] = 2*ones(51,N_RPS[i],1,10)-2*load(path_RPS*f,"data")
		end
	end
end

# ╔═╡ 81686c37-e502-4172-a163-8fbf43791382
md"## |11....11> normalized"

# ╔═╡ 854c763a-4645-4b02-985f-10e23f461350
begin
	path_PSI0 = pwd()*"/PSI0/"
	F_PSI0 = 13
	jobids_PSI0 = Vector{String}(undef,F_PSI0)
	params_PSI0 = Vector{Any}(undef,F_PSI0)
	data_PSI0 = Vector{Array{Float64,4}}(undef,F_PSI0)
	N_PSI0 = [8,9,10,11,12,13,14,15,16,17,18,19,20]
	files_PSI0 = ["5203852_N8.jld2","5203853_N9.jld2","5203854_N10.jld2","5203855_N11.jld2","5203857_N12.jld2","5203858_N13.jld2","5203859_N14.jld2","5203860_N15.jld2","5203861_N16.jld2","5203862_N17.jld2","5203863_N18.jld2","5203864_N19.jld2","5203856_N20.jld2"]
	for (i,f) in enumerate(files_PSI0)
		jobids_PSI0[i] = load(path_PSI0*f,"jobid")
		params_PSI0[i] = load(path_PSI0*f,"params")
		data_PSI0[i] = 2*ones(51,N_PSI0[i],1,1)-2*load(path_PSI0*f,"data")
	end
end

# ╔═╡ 28491406-99a1-11ec-3ca4-ef0928b5ed59
begin
	function state_mean(A,n_states=10)
		return mean(A[:,:,:,1:n_states];dims=4)[:,:,:,1]
	end
	
	function state_std(A,n_states)
		return std(A[:,:,:,1:n_states];dims=4)[:,:,:,1]
	end

	function pos_mean(A)
		return mean(A;dims=2)[:,1]
	end
	
	function reduce_by_last(A)
		return A[:,:,1]
	end
end

# ╔═╡ 65b15201-1b16-49ba-a3ed-7416bf43e206
md"# Qualitative Light Cones"

# ╔═╡ 4b448258-8848-43d0-a1ca-d94be49f69e4
@bind idx_RS Slider(1:F_RS)

# ╔═╡ d5a0df00-c519-4f43-ad79-5b24405726a0
@bind idx_PSI0 Slider(1:F_PSI0)

# ╔═╡ 98f4b267-2171-4596-bc14-9c150d0bae74
@bind idx_RPS Slider(1:F_RPS)

# ╔═╡ 8c447adc-8212-4cba-9394-3a2eb03a8d79
@bind si Slider(1:10)

# ╔═╡ 4ad29cc2-4558-4618-9e7f-585737522cd9
begin
	plot_RS = heatmap(1:N_RS[idx_RS],0:0.1:5,data_RS[idx_RS][:,:,1,si];title="RS for N=$(N_RS[idx_RS])",c=:viridis)
	
	plot_RPS = heatmap(1:N_RPS[idx_RPS],0:0.1:5,data_RPS[idx_RPS][:,:,1,si];title="RPS for N=$(N_RPS[idx_RPS])",c=:viridis)

	plot_PSI0 = heatmap(1:N_PSI0[idx_PSI0],0:0.1:5,data_PSI0[idx_PSI0][:,:,1,1];title="PSI0 for N=$(N_PSI0[idx_PSI0])",c=:viridis)
	
	
	plot(plot_RS, plot_RPS, plot_PSI0, layout = (1, 3), legend = false)
	#latexstring(raw"<|[σ_i,σ_j]|^2>\mathrm{\;for\;}"*"N={$(N[idx])}") #Note: Raise issue as this does not work in title!
end

# ╔═╡ 6bd295fd-269c-45f2-b0d5-50640fd03bbb
@bind idx_L Slider(1:F_L)

# ╔═╡ 8d2670ab-8e09-451c-9dc3-f07b39ded685
@bind si_L Slider(1:100)

# ╔═╡ cd527e24-c144-4ffc-9584-d5d9a4f0c0a0
begin
	plot_RS_L = heatmap(1:N_L[idx_L],0:0.1:5,data_L[idx_L][:,:,1,si_L];title="RS_L for N=$(N_L[idx_L])",c=:viridis)
	
	plot_RPS_L = heatmap(1:N_RPS[idx_RPS],0:0.1:5,data_RPS[idx_RPS][:,:,1,si_L];title="RPS for N=$(N_RPS[idx_RPS])",c=:viridis)

	plot(plot_RS_L, plot_RPS_L, layout = (1, 2), legend = false)
	#latexstring(raw"<|[σ_i,σ_j]|^2>\mathrm{\;for\;}"*"N={$(N[idx])}") #Note: Raise issue as this does not work in title!
end

# ╔═╡ 84cbed96-bf93-4650-a053-8fa950635c8e
@bind n_states Slider(2:10)

# ╔═╡ 6c52bd5c-a2c2-4013-b0e0-cfcffdc8cd22
@bind n_states_l Slider(2:100)

# ╔═╡ 0e75b8de-0df7-4a7b-b0dc-214486c3f08e
begin
	data_mean_RS = Vector{Array{Float64,2}}(undef,F_RS)
	data_std_RS = Vector{Array{Float64,2}}(undef,F_RS)
	for i in 1:F_RS
		data_mean_RS[i] = reduce_by_last(state_mean(data_RS[i],n_states))
		data_std_RS[i] = reduce_by_last(state_std(data_RS[i],n_states))
	end

	data_mean_RPS = Vector{Array{Float64,2}}(undef,F_RPS)
	data_std_RPS = Vector{Array{Float64,2}}(undef,F_RPS)
	for i in 1:F_RPS
		data_mean_RPS[i] = reduce_by_last(state_mean(data_RPS[i],n_states_l))
		data_std_RPS[i] = reduce_by_last(state_std(data_RPS[i],n_states_l))
	end

	data_mean_L = Vector{Array{Float64,2}}(undef,F_L)
	data_std_L = Vector{Array{Float64,2}}(undef,F_L)
	for i in 1:F_L
		data_mean_L[i] = reduce_by_last(state_mean(data_L[i],n_states_l))
		data_std_L[i] = reduce_by_last(state_std(data_L[i],n_states_l))
	end

	data_mean_PSI0 = Vector{Array{Float64,2}}(undef,F_PSI0)
	data_std_PSI0 = Vector{Array{Float64,2}}(undef,F_PSI0)
	for i in 1:F_PSI0
		data_mean_PSI0[i] = reduce_by_last(state_mean(data_PSI0[i],1))
		data_std_PSI0[i] = reduce_by_last(state_std(data_PSI0[i],1))
	end
end

# ╔═╡ 0ba548f8-941f-4173-a58b-6e3348e44208
begin
	plot_RS_mean = heatmap(1:N_RS[idx_RS],0:0.1:5,data_mean_RS[idx_RS][:,:];title="RS for N=$(N_RS[idx_RS])",c=:viridis)
	
	plot_RPS_mean = heatmap(1:N_RPS[idx_RPS],0:0.1:5,data_mean_RPS[idx_RPS][:,:];title="RPS for N=$(N_RPS[idx_RPS])",c=:viridis)

	plot_PSI0_mean = heatmap(1:N_PSI0[idx_PSI0],0:0.1:5,data_mean_PSI0[idx_PSI0][:,:];title="PSI0 for N=$(N_PSI0[idx_PSI0])",c=:viridis)
	
	
	plot(plot_RS_mean, plot_RPS_mean, plot_PSI0_mean, layout = (1, 3), legend = false)
	#latexstring(raw"<|[σ_i,σ_j]|^2>\mathrm{\;for\;}"*"N={$(N[idx])}") #Note: Raise issue as this does not work in title!
end

# ╔═╡ 9ad60e8b-e180-4126-ad79-2a908ee61cc9
begin
	plot_RS_L_mean = heatmap(1:N_L[idx_L],0:0.1:5,data_mean_L[idx_L][:,:];title="RS_L for N=$(N_L[idx_L])",c=:viridis)
	
	plot_RPS_L_mean = heatmap(1:N_RPS[idx_RPS],0:0.1:5,data_mean_RPS[idx_RPS][:,:];title="RPS for N=$(N_RPS[idx_RPS])",c=:viridis)

	plot(plot_RS_L_mean, plot_RPS_L_mean, layout = (1, 2), legend = false)
	#latexstring(raw"<|[σ_i,σ_j]|^2>\mathrm{\;for\;}"*"N={$(N[idx])}") #Note: Raise issue as this does not work in title!
end

# ╔═╡ da45af9d-c845-42ae-b9d2-a5a02efd5a54
md"# OTOCs at specific positions"

# ╔═╡ 2c045bae-a1dc-4644-842a-70f7c762655a
@bind idx_PSI02 Slider(1:F_PSI0)

# ╔═╡ 18942f53-a04b-4a6c-9db8-5869a0e129f4
@bind idx_RPS2 Slider(1:F_RPS)

# ╔═╡ ec7467d9-7823-46c4-84aa-c02934c03f38
@bind idx_RS2 Slider(1:F_RS)

# ╔═╡ 171e7121-cc30-477d-922b-b6227edbb09b
N_RS[idx_RS2]

# ╔═╡ 0e2a62cc-91e3-4730-b6fc-6baf7efd4aaa
@bind pos Slider(1:18)

# ╔═╡ dc9b1db8-bf3d-4e62-b9cc-eb0b37f84eb7
begin
	plot(0:0.1:5,data_mean_RS[idx_RS2][1:51,pos],yerrors=data_std_RS[idx_RS2][1:51,pos]./sqrt(10),title="Comparison for N=20, i=3, j=$(pos)",xlabel="t",ylabel="<|[σ_i(t),σ_j]|²>",label="10 Haar Random States: N=$(N_RS[idx_RS2]),i=3,j=$(pos)")
	
	plot!(0:0.1:5,data_mean_RPS[idx_RPS2][1:51,pos],ribbon=data_std_RPS[idx_RPS2][1:51,pos]./sqrt(n_states_l),xlabel="t",ylabel="<|[σ_i(t),σ_j]|²>",label="100 Random Product States: N=$(N_RPS[idx_RPS2]),i=3,j=$(pos)")

	plot!(0:0.1:5,data_mean_PSI0[idx_PSI02][1:51,pos],yerrors=data_std_PSI0[idx_PSI02][1:51,pos],xlabel="t",ylabel="<|[σ_i(t),σ_j]|²>",label="|Ψ0> = |11..11> N=$(N_PSI0[idx_PSI02]),i=3,j=$(pos)",legend=:bottomright)
	
	#ribbon for shaded area
end

# ╔═╡ 2082a75b-217d-44d0-a5da-40327434df67
begin
	plot(0.1:0.1:5,data_mean_RS[idx_RS2][2:51,pos],yerrors=data_std_RS[idx_RS2][2:51,pos]./sqrt(10),title="Comparison for N<=20 (100 states)",xlabel="t",ylabel="<|[σ_i(t),σ_j]|²>",label="RS: N=$(N_RS[idx_RS2]),i=3,j=$(pos)",yaxis=:log,xaxis=:log,legend=:bottomright)
	
	plot!(0.1:0.1:5,data_mean_RPS[idx_RPS2][2:51,pos],ribbon=data_std_RPS[idx_RPS2][2:51,pos]./sqrt(100),xlabel="t",ylabel="<|[σ_i(t),σ_j]|²>",label="RPS: N=$(N_RPS[idx_RPS2]),i=3,j=$(pos)",yaxis=:log,xaxis=:log)

	plot!(0.1:0.1:5,data_mean_PSI0[idx_PSI02][2:51,pos],yerrors=data_std_PSI0[idx_PSI02][2:51,pos],xlabel="t",ylabel="<|[σ_i(t),σ_j]|²>",label="Ψ0:N=$(N_PSI0[idx_PSI02]),i=3,j=$(pos)",yaxis=:log,xaxis=:log)
	
	#ribbon for shaded area
end

# ╔═╡ 27860efe-48f0-4547-800f-a6a5daf105a4
#Choose files with smaller RS samples for this comparison or redo as below

# ╔═╡ 74dabc7f-d6a1-43d0-8062-6d7fbb6c1e45
@bind idx_RS4 Slider(1:F_RS)

# ╔═╡ 19a3c921-576d-4e1c-9147-d14057de522f
@bind idx_L4 Slider(1:F_L)

# ╔═╡ 2e24dd9f-632f-4307-914d-c42082f72e68
@bind pos4 Slider(1:18)

# ╔═╡ 561dd738-381c-439f-a191-4d371f566597
begin
	plot(0:0.1:5,data_mean_RS[idx_RS4][1:51,pos4],yerrors=data_std_RS[idx_RS4][1:51,pos4],title="Random States, varying sample sizes (10 vs 100)",xlabel="t",ylabel="<|[σ_i(t),σ_j]|²>",label="RS N=$(N_RS[idx_RS4]),i=3,j=$(pos4)")
	
	plot!(0:0.1:5,data_mean_L[idx_L4][1:51,pos4],ribbon=data_std_L[idx_L4][1:51,pos4],xlabel="t",ylabel="<|[σ_i(t),σ_j]|²>",label="RS L: N=$(N_L[idx_L4]),i=3,j=$(pos4)")
end

# ╔═╡ 71ee2089-eb23-4deb-8682-e7c11a94a083
begin
	plot(0.1:0.1:5,data_mean_RS[idx_RS4][2:51,pos4],yerrors=data_std_RS[idx_RS4][2:51,pos4],title="Random States, varying sample sizes (10 vs 100)",xlabel="t",ylabel="<|[σ_i(t),σ_j]|²>",label="RS N=$(N_RS[idx_RS4]),i=3,j=$(pos4)",yaxis=:log,xaxis=:log)
	
	plot!(0.1:0.1:5,data_mean_L[idx_L4][2:51,pos4],ribbon=data_std_L[idx_L4][2:51,pos4],xlabel="t",ylabel="<|[σ_i(t),σ_j]|²>",label="RS L: N=$(N_L[idx_L4]),i=3,j=$(pos4)",yaxis=:log,xaxis=:log)
end

# ╔═╡ 8bf67ebf-c341-4c86-880f-aa30fd1e7909
#Now different analysis

# ╔═╡ c77606b9-3bef-4cc7-8e13-f6d3951e0746
@bind n_states_RPS_s Slider(1:10)

# ╔═╡ 0ac8c8a3-488a-4e15-ae77-1e1c1dd75ad0
@bind n_states_RPS_l Slider(1:100)

# ╔═╡ aa13c33c-e214-4e8f-a7c4-90a6ab5129d7
@bind idx_RPS5 Slider(1:13)

# ╔═╡ f0ae6427-b84b-4337-815f-58562dff26b6
@bind idx_RS5 Slider(1:13)

# ╔═╡ 36353cad-4ee8-4576-ab02-bf0c650bb0f9
@bind pos5 Slider(1:18)

# ╔═╡ 04bc256f-3e8b-4878-b937-1ad4d058e189
begin
	data_mean_RPS5 = Vector{Array{Float64,2}}(undef,13)
	data_std_RPS5 = Vector{Array{Float64,2}}(undef,13)
	data_mean_RPS6 = Vector{Array{Float64,2}}(undef,13)
	data_std_RPS6 = Vector{Array{Float64,2}}(undef,13)
	for i in 1:13
		data_mean_RPS5[i] = reduce_by_last(state_mean(data_RPS[i],n_states_RPS_s))
		data_std_RPS5[i] = reduce_by_last(state_std(data_RPS[i],n_states_RPS_s))
		data_mean_RPS6[i] = reduce_by_last(state_mean(data_RPS[i],n_states_RPS_l))
		data_std_RPS6[i] = reduce_by_last(state_std(data_RPS[i],n_states_RPS_l))
	end
end

# ╔═╡ d0a9fb36-ca76-4c1f-8201-11a88cee074d
begin
	plot(0:0.1:5,data_mean_RPS5[idx_RPS5][1:51,pos5],yerrors=data_std_RPS5[idx_RPS5][1:51,pos5]./sqrt(n_states_RPS_s),title="Random Product States, var. sample sizes ($(n_states_RPS_s) vs. $(n_states_RPS_l))",xlabel="t",ylabel="<|[σ_i(t),σ_j]|²>",label="RPS 10: N=$(N_RPS[idx_RPS5]),i=3,j=$(pos5)")
	
	plot!(0:0.1:5,data_mean_RPS6[idx_RPS5][1:51,pos5],ribbon=data_std_RPS6[idx_RPS5][1:51,pos5]./sqrt(n_states_RPS_l),xlabel="t",ylabel="<|[σ_i(t),σ_j]|²>",label="RPS 100: N=$(N_RPS[idx_RPS5]),i=3,j=$(pos5)")

	plot!(0:0.1:5,data_mean_RS[idx_RS5][1:51,pos5],ribbon=data_std_RS[idx_RS5][1:51,pos5]./sqrt(10),xlabel="t",ylabel="<|[σ_i(t),σ_j]|²>",label="RS N=$(N_RS[idx_RS5]),i=3,j=$(pos5)")
end

# ╔═╡ 90a64aa8-3901-4f20-aae7-eb984a01a13a
begin
	plot(0.1:0.1:5,data_mean_RPS5[idx_RPS5][2:51,pos5],yerrors=data_std_RPS5[idx_RPS5][2:51,pos5]./sqrt(n_states_RPS_s),title="Random Product States, var. sample sizes ($(n_states_RPS_s) vs. $(n_states_RPS_l))",xlabel="t",ylabel="<|[σ_i(t),σ_j]|²>",label="RPS N=$(N_RPS[idx_RPS5]),i=3,j=$(pos5)",yaxis=:log,xaxis=:log)
	
	plot!(0.1:0.1:5,data_mean_RPS6[idx_RPS5][2:51,pos5],ribbon=data_std_RPS6[idx_RPS5][2:51,pos5]./sqrt(n_states_RPS_l),xlabel="t",ylabel="<|[σ_i(t),σ_j]|²>",label="RPS L: N=$(N_RPS[idx_RPS5]),i=3,j=$(pos5)",yaxis=:log,xaxis=:log)

	plot!(0.1:0.1:5,data_mean_RS[idx_RS5][2:51,pos5],ribbon=data_std_RS[idx_RS5][2:51,pos5]./sqrt(10),xlabel="t",ylabel="<|[σ_i(t),σ_j]|²>",label="RS N=$(N_RS[idx_RS5]),i=3,j=$(pos5)",yaxis=:log,xaxis=:log,legend=true)
end

# ╔═╡ 26218b0a-035d-477f-8e27-7f5a11d80771
md"# Scaling of convergence to mean"

# ╔═╡ 55526188-6e02-4055-9575-fd640d0e17ac
@bind n_states_RPS_c Slider(1:100)

# ╔═╡ a4e13f59-3446-405d-86c1-995c82fa3593
@bind idx_RPS7 Slider(1:13)

# ╔═╡ 1e3e9f8a-d603-484b-a8cf-78a33fc2773a
N_RPS[idx_RPS7]

# ╔═╡ a95bfa09-05fc-4da6-98bf-5c7d4f6f1cd1
@bind idx_RS7 Slider(1:13)

# ╔═╡ b2e5934a-090c-4110-95cb-b571d237a363
N_RS[idx_RS7]

# ╔═╡ dd19bdc1-8d1d-4f20-9f1c-d9a98e2f040e


# ╔═╡ c1398f9f-3990-492a-b21e-c9d8e18fdd3c
@bind pos7 Slider(1:17)

# ╔═╡ 1fcb7174-f429-455a-855c-a82d0344fe46
begin
	data_mean_RPS7 = Vector{Array{Float64,2}}(undef,13)
	data_std_RPS7 = Vector{Array{Float64,2}}(undef,13)
	data_mean_RPS8 = Vector{Array{Float64,2}}(undef,13)
	data_std_RPS8 = Vector{Array{Float64,2}}(undef,13)
	for i in 1:13
		data_mean_RPS7[i] = reduce_by_last(state_mean(data_RPS[i],n_states_RPS_c))
		data_std_RPS7[i] = reduce_by_last(state_std(data_RPS[i],n_states_RPS_c))
		data_mean_RPS8[i] = reduce_by_last(state_mean(data_RPS[i],100))
		data_std_RPS8[i] = reduce_by_last(state_std(data_RPS[i],100))
	end
end

# ╔═╡ 0189c9ad-6891-4152-84b5-2bfabfab69c7
begin
	σμ = zeros(13,20)
	N = [8,9,10,11,12,13,14,15,16,17,18,19,20]
	for i in 1:13
		for p in 1:N[i]
			σμ[i,p] = sum(data_std_RPS7[i][:,p])./sqrt(n_states_RPS_c)/51
		end
	end
end

# ╔═╡ 60804ae2-2fb3-4a85-9a24-91cf1d7777f3
begin
	plot(0:0.1:5,data_std_RPS7[1][:,pos7]./sqrt(n_states_RPS_c),label="N=$(N[1])",title="Standardfehler Mittewert i=3,j=$(pos7)",xlabel="t")
	for i in 2:9
		plot!(0:0.1:5,data_std_RPS7[i][:,pos7]./sqrt(n_states_RPS_c),label="N=$(N[i])")
	end
	plot!(0:0.1:5,data_std_RPS7[10][:,pos7]./sqrt(n_states_RPS_c),label="N=$(N[10])")

end

# ╔═╡ 175816a8-b892-4f41-b411-63612c80c8fe
plot(N,σμ,ylims=(-0.1,0.4),legend=false)

# ╔═╡ 3adb899a-5168-4209-a8c4-43687e61d47c
md"# Deviation from RS mean"

# ╔═╡ f971aa47-83d7-4d5c-8edc-5e2a9c95594c
begin
	plot(0:0.1:5,data_mean_RPS7[idx_RPS7][1:51,pos7]-data_mean_RS[idx_RS7][1:51,pos7],yerrors=data_std_RPS7[idx_RPS7][1:51,pos7]./sqrt(n_states_RPS_c),title="Deviation from RS mean (#Shots = $(n_states_RPS_c))",xlabel="t",ylabel="m-μ",label="RPS-RS N=$(N_RPS[idx_RPS7]),i=3,j=$(pos7)",ylims=(-0.1,0.1))
end

# ╔═╡ 9aa1b1d6-06e3-400e-bc59-aa6a113196ac
begin
	plot(0.1:0.1:5,(data_mean_RPS7[idx_RPS7][2:51,pos7]-data_mean_RS[idx_RS7][2:51,pos7])./data_mean_RS[idx_RS7][2:51,pos7]*100,yerrors=data_std_RPS7[idx_RPS7][2:51,pos7]./sqrt(n_states_RPS_c)./data_mean_RS[idx_RS7][2:51,pos7]*100,title="Deviation from RS mean (#Shots = $(n_states_RPS_c))",xlabel="t",ylabel="(m-μ)/μ [%]",label="RPS-RS N=$(N_RPS[idx_RPS7]),i=3,j=$(pos7)",ylims=(-10,10))
end

# ╔═╡ b4a557ce-62dc-41d0-a0f6-b76ff35212aa
#Student t-test

# ╔═╡ 7c8d2139-5ed3-4f88-b439-f019ae94b423


# ╔═╡ a904d996-beb6-43e4-900a-a62f46566d2a
plot(0:0.1:5,data_std[idx2][:,pos],title="Absolute typicality std for N=$(N[idx2]), n_states=$(n_states)",label="i=3,j=$(pos)",xlabel="t",ylabel="std") #ribbon for

# ╔═╡ 75532f27-8989-4ed7-b206-a45135b5eecb
plot(0.1:0.1:5,data_std[idx2][2:51,pos]./data_mean[idx2][2:51,pos]*100,title="Relative typicality std for N=$(N[idx2]), n_states=$(n_states)",label="i=3,j=$(pos)",xlabel="t",ylabel="std / %")

# ╔═╡ 1eaef926-fecd-4a37-b5db-9c706e2e16d7
plot(0:0.1:5,pos_mean(data_std[idx2]),title="Absolute typicality std for N=$(N[idx2]), n_states=$(n_states)",label="i=3",xlabel="t",ylabel="std") #ribbon for

# ╔═╡ 6e87bf6f-c7b0-4f67-ae84-820f3e40228b
plot(0.1:0.1:5,pos_mean(data_std[idx2][2:51,:])./pos_mean(data_mean[idx2][2:51,:])*100,title="Relative typicality std for N=$(N[idx2]), n_states=$(n_states)",label="i=3",xlabel="t",ylabel="std / %")

# ╔═╡ abf982fb-2a9b-41f9-b31a-be3921429efc
0

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
JLD2 = "~0.4.22"
Plots = "~1.26.0"
PlutoUI = "~0.7.36"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c9a6160317d1abe9c44b3beb367fd448117679ca"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.13.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "12fc73e5e0af68ad3137b886e3f7c1eacfca2640"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.17.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "96b0bc6c52df76506efc8a441c6cf1adcb1babc4"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.42.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ae13fcbc7ab8f16b0856729b050ef0c446aa3492"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.4+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "80ced645013a5dbdc52cf70329399c35ce007fae"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.13.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "9f836fb62492f4b0f0d3b06f55983f2704ed0883"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.0"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a6c850d77ad5118ad3be4bd188919ce97fffac47"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.64.0+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "81b9477b49402b47fbe7f7ae0b252077f53e4a08"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.22"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a6552bfeab40de157a297d84e03ade4b8177677f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.12"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "3f7cb7157ef860c637f3f4929c8ed5d9716933c6"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.7"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
git-tree-sha1 = "737a5957f387b17e74d4ad2f440eb330b39a62c5"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "648107615c15d4e09f7eca16307bc821c1f718d8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.13+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "85b5da0fa43588c75bb1ff986493443f821c70b7"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "6f1b25e8ea06279b5689263cc538f51331d7ca17"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.1.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "23d109aad5d225e945c813c6ebef79104beda955"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.26.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "2c87c85e397b7ffed5ffec054f532d4edd05d901"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.36"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "de893592a221142f3db370f48290e3a2ef39998f"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.4"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "c6c0f690d0cc7caddb74cef7aa847b824a16b256"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "995a812c6f7edea7527bb570f0ac39d0fb15663c"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.1"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "74fb527333e72ada2dd9ef77d98e4991fb185f04"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.1"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c3d8ba7f3fa0625b062b82853a7d5229cb728b6b"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.1"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "57617b34fa34f91d536eb265df67c2d4519b8b98"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.5"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╠═c39afa38-c3d8-49cd-b730-3919e7b7a069
# ╠═18a0434c-bacf-411f-a840-fb69a9a91c7c
# ╠═6aec44d6-1e21-4f5d-85f8-a25cd942e0c4
# ╠═6612a15a-c3ca-49da-823f-808e8fc6f8ac
# ╠═2e9d9a87-b8ee-48f4-bfcd-33514a9c0250
# ╠═6f6f4f0f-d83d-48c8-afed-9b2ba5c02523
# ╠═e662f523-7129-4ea7-8564-ab2473449b7c
# ╠═a52bd686-43c8-49af-993d-174a5d10119b
# ╠═79d282be-08d3-4bf3-b700-eee967e1815a
# ╠═f636cf34-b778-4259-abd1-bede4cf91a0f
# ╠═64ac6e06-1c34-4430-9bc4-d0bb9097834f
# ╠═397f89e9-2051-4109-90ac-47e31a993637
# ╠═81686c37-e502-4172-a163-8fbf43791382
# ╠═854c763a-4645-4b02-985f-10e23f461350
# ╠═28491406-99a1-11ec-3ca4-ef0928b5ed59
# ╠═65b15201-1b16-49ba-a3ed-7416bf43e206
# ╠═4b448258-8848-43d0-a1ca-d94be49f69e4
# ╠═d5a0df00-c519-4f43-ad79-5b24405726a0
# ╠═98f4b267-2171-4596-bc14-9c150d0bae74
# ╠═8c447adc-8212-4cba-9394-3a2eb03a8d79
# ╠═4ad29cc2-4558-4618-9e7f-585737522cd9
# ╠═6bd295fd-269c-45f2-b0d5-50640fd03bbb
# ╠═8d2670ab-8e09-451c-9dc3-f07b39ded685
# ╠═cd527e24-c144-4ffc-9584-d5d9a4f0c0a0
# ╠═84cbed96-bf93-4650-a053-8fa950635c8e
# ╠═6c52bd5c-a2c2-4013-b0e0-cfcffdc8cd22
# ╠═0e75b8de-0df7-4a7b-b0dc-214486c3f08e
# ╠═0ba548f8-941f-4173-a58b-6e3348e44208
# ╠═9ad60e8b-e180-4126-ad79-2a908ee61cc9
# ╠═da45af9d-c845-42ae-b9d2-a5a02efd5a54
# ╠═2c045bae-a1dc-4644-842a-70f7c762655a
# ╠═18942f53-a04b-4a6c-9db8-5869a0e129f4
# ╠═ec7467d9-7823-46c4-84aa-c02934c03f38
# ╠═171e7121-cc30-477d-922b-b6227edbb09b
# ╠═0e2a62cc-91e3-4730-b6fc-6baf7efd4aaa
# ╠═dc9b1db8-bf3d-4e62-b9cc-eb0b37f84eb7
# ╠═2082a75b-217d-44d0-a5da-40327434df67
# ╠═27860efe-48f0-4547-800f-a6a5daf105a4
# ╠═74dabc7f-d6a1-43d0-8062-6d7fbb6c1e45
# ╠═19a3c921-576d-4e1c-9147-d14057de522f
# ╠═2e24dd9f-632f-4307-914d-c42082f72e68
# ╠═561dd738-381c-439f-a191-4d371f566597
# ╠═71ee2089-eb23-4deb-8682-e7c11a94a083
# ╠═8bf67ebf-c341-4c86-880f-aa30fd1e7909
# ╠═c77606b9-3bef-4cc7-8e13-f6d3951e0746
# ╠═0ac8c8a3-488a-4e15-ae77-1e1c1dd75ad0
# ╠═aa13c33c-e214-4e8f-a7c4-90a6ab5129d7
# ╠═f0ae6427-b84b-4337-815f-58562dff26b6
# ╠═36353cad-4ee8-4576-ab02-bf0c650bb0f9
# ╠═04bc256f-3e8b-4878-b937-1ad4d058e189
# ╠═d0a9fb36-ca76-4c1f-8201-11a88cee074d
# ╠═90a64aa8-3901-4f20-aae7-eb984a01a13a
# ╠═26218b0a-035d-477f-8e27-7f5a11d80771
# ╠═55526188-6e02-4055-9575-fd640d0e17ac
# ╠═a4e13f59-3446-405d-86c1-995c82fa3593
# ╠═1e3e9f8a-d603-484b-a8cf-78a33fc2773a
# ╠═a95bfa09-05fc-4da6-98bf-5c7d4f6f1cd1
# ╠═b2e5934a-090c-4110-95cb-b571d237a363
# ╠═dd19bdc1-8d1d-4f20-9f1c-d9a98e2f040e
# ╠═c1398f9f-3990-492a-b21e-c9d8e18fdd3c
# ╠═1fcb7174-f429-455a-855c-a82d0344fe46
# ╠═0189c9ad-6891-4152-84b5-2bfabfab69c7
# ╠═60804ae2-2fb3-4a85-9a24-91cf1d7777f3
# ╠═175816a8-b892-4f41-b411-63612c80c8fe
# ╠═3adb899a-5168-4209-a8c4-43687e61d47c
# ╠═f971aa47-83d7-4d5c-8edc-5e2a9c95594c
# ╠═9aa1b1d6-06e3-400e-bc59-aa6a113196ac
# ╠═b4a557ce-62dc-41d0-a0f6-b76ff35212aa
# ╠═7c8d2139-5ed3-4f88-b439-f019ae94b423
# ╠═a904d996-beb6-43e4-900a-a62f46566d2a
# ╠═75532f27-8989-4ed7-b206-a45135b5eecb
# ╠═1eaef926-fecd-4a37-b5db-9c706e2e16d7
# ╠═6e87bf6f-c7b0-4f67-ae84-820f3e40228b
# ╠═abf982fb-2a9b-41f9-b31a-be3921429efc
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
