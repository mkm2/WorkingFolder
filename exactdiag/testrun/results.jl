### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ 51dbff79-d4aa-4a47-920b-9c481d5c25fa
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 514bb0f4-3b25-11ed-1979-71d8df805c28
using LinearAlgebra,Plots,JLD2,Statistics,PlutoUI

# ╔═╡ 559a03c0-a43a-4fa5-a558-5b1ae71dd7f9
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

# ╔═╡ 757cf2b6-3fb0-4f8a-9977-111059553d13
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

# ╔═╡ 2cea701f-f745-4d29-bcba-5983e10182ba
trange = logrange(-2,4,1e4)

# ╔═╡ a812e486-b2d1-4172-890c-b92524b639d1
length(trange)

# ╔═╡ c7201d3b-ac4b-45ff-98ac-606e23003e5a
md"# h = 0"

# ╔═╡ b3373faa-a84a-404c-8f38-c89e28fad71b
N = 15

# ╔═╡ 6436b188-b3f4-4f65-aef9-6990c2613f5d
begin
	path_0 = pwd()*"/h=0/"
	#F_RS = 10
	f_0 = "7125541_N15_ED.jld2"
	#combine_files(["6192397_N19_RS.jld2","6192398_N19_RS.jld2","6192399_N19_RS.jld2","6192400_N19_RS.jld2","6192401_N19_RS.jld2","6192402_N19_RS.jld2","6192403_N19_RS.jld2","6192404_N19_RS.jld2","6192405_N19_RS.jld2","6192400_N19_RS.jld2"],path_RS,f_RS)
	
	jobids_0 = load(path_0*f_0,"jobid")
	params_0 = load(path_0*f_0,"params")
	data_0 = 2*ones(56,N,50)-2*load(path_0*f_0,"data")
	
	size(data_0)
end

# ╔═╡ e255cbdf-6794-460c-8b72-98e1e6f68fbf
begin
	data_0_mean = disorder_mean(data_0,size(data_0)[3])
	data_0_std = disorder_std(data_0,size(data_0)[3])
end

# ╔═╡ 3449b6fa-ba13-414c-8823-287d69c8adf7
plot(trange[2:56],data_0_mean[2:56,:],xaxis=:log,legend=nothing,ribbon=data_0_std[2:56,:]/sqrt(50))

# ╔═╡ 5ee2f089-db1a-471b-b3e5-c259053cbbd2
heatmap(1:15,trange[2:56],data_0_mean[2:56,:],yaxis=:log)

# ╔═╡ ff4eb17b-bf50-42b6-ac95-69a681e35934
md"# h=4"

# ╔═╡ 459fb66e-ba82-4cd6-88b9-19ca3d395258
shots4 = 1000

# ╔═╡ e5682b15-3d82-4734-9f1f-100b0c4baa58
begin
	path_4 = pwd()*"/h=4/"
	f_4 = "7125556-7127856_N15_ED.jld2"
	combine_files(["7125556_N15_ED.jld2","7125557_N15_ED.jld2","7125558_N15_ED.jld2","7125559_N15_ED.jld2","7125560_N15_ED.jld2","7125561_N15_ED.jld2","7125562_N15_ED.jld2","7125564_N15_ED.jld2","7125565_N15_ED.jld2","7125566_N15_ED.jld2","7125567_N15_ED.jld2","7125568_N15_ED.jld2","7125569_N15_ED.jld2","7125570_N15_ED.jld2","7125571_N15_ED.jld2","7125572_N15_ED.jld2","7125573_N15_ED.jld2","7125574_N15_ED.jld2","7125575_N15_ED.jld2","7127856_N15_ED.jld2"],path_4,f_4)
	
	jobids_4 = load(path_4*f_4,"jobid")
	params_4 = load(path_4*f_4,"params")
	data_4 = 2*ones(56,N,1000)-2*cat(load(path_4*f_4,"data")...,dims=3)

	size(data_4)
end

# ╔═╡ 4a1b1d89-c810-416b-a361-a33fc5542200
begin
	data_4_mean = disorder_mean(data_4,size(data_4)[3])
	data_4_std = disorder_std(data_4,size(data_4)[3])
end

# ╔═╡ 61baa510-4196-463f-bcac-dafd78785967
plot(trange[2:56],data_4_mean[2:56,:],xaxis=:log,legend=nothing,ribbon=data_4_std[2:56,:]/sqrt(shots4))

# ╔═╡ 45b4700a-399a-42c2-ac5a-07fdbe80195d
heatmap(1:15,trange[2:56],data_4_mean[2:56,:],yaxis=:log)

# ╔═╡ f8476e50-cdc7-4239-ace4-cfc90fd94e4e
md"# h=12"

# ╔═╡ 033e3abe-318b-4ccc-a18d-f4b692312d92
shots12 = 100

# ╔═╡ d627fd5a-cddd-4179-8796-4bc967585ead
begin
	path_12 = pwd()*"/h=12/"
	f_12 = "7125576-7125595_N15_ED.jld2"
	combine_files(["7125576_N15_ED.jld2","7125577_N15_ED.jld2","7125578_N15_ED.jld2","7125579_N15_ED.jld2","7125580_N15_ED.jld2","7125581_N15_ED.jld2","7125582_N15_ED.jld2","7125583_N15_ED.jld2","7125584_N15_ED.jld2","7125585_N15_ED.jld2","7125586_N15_ED.jld2","7125587_N15_ED.jld2","7125588_N15_ED.jld2","7125589_N15_ED.jld2","7125590_N15_ED.jld2","7125591_N15_ED.jld2","7125592_N15_ED.jld2","7125593_N15_ED.jld2","7125594_N15_ED.jld2","7125595_N15_ED.jld2"],path_12,f_12)
	
	jobids_12 = load(path_12*f_12,"jobid")
	params_12 = load(path_12*f_12,"params")
	data_12 = 2*ones(56,N,1000)-2*cat(load(path_12*f_12,"data")...,dims=3)

	size(data_12)
end

# ╔═╡ 42de57f4-9742-4e13-b28c-ed9b6aa6d8be
begin
	data_12_mean = disorder_mean(data_12,shots12)
	data_12_std = disorder_std(data_12,shots12)
end

# ╔═╡ b13f1b95-68ee-4b24-8584-e7a837ab0b5d
plot(trange[2:56],data_12_mean[2:56,:],xaxis=:log,legend=nothing,ribbon=data_12_std[2:56,:]/sqrt(shots12))

# ╔═╡ 0c26b8ce-6ac8-401e-951f-a1fe9f6d272a
heatmap(1:15,trange[2:56],data_12_mean[2:56,:],yaxis=:log)

# ╔═╡ Cell order:
# ╠═514bb0f4-3b25-11ed-1979-71d8df805c28
# ╠═51dbff79-d4aa-4a47-920b-9c481d5c25fa
# ╠═559a03c0-a43a-4fa5-a558-5b1ae71dd7f9
# ╠═757cf2b6-3fb0-4f8a-9977-111059553d13
# ╠═2cea701f-f745-4d29-bcba-5983e10182ba
# ╠═a812e486-b2d1-4172-890c-b92524b639d1
# ╠═c7201d3b-ac4b-45ff-98ac-606e23003e5a
# ╠═b3373faa-a84a-404c-8f38-c89e28fad71b
# ╠═6436b188-b3f4-4f65-aef9-6990c2613f5d
# ╠═e255cbdf-6794-460c-8b72-98e1e6f68fbf
# ╠═3449b6fa-ba13-414c-8823-287d69c8adf7
# ╠═5ee2f089-db1a-471b-b3e5-c259053cbbd2
# ╠═ff4eb17b-bf50-42b6-ac95-69a681e35934
# ╠═459fb66e-ba82-4cd6-88b9-19ca3d395258
# ╠═e5682b15-3d82-4734-9f1f-100b0c4baa58
# ╠═4a1b1d89-c810-416b-a361-a33fc5542200
# ╠═61baa510-4196-463f-bcac-dafd78785967
# ╠═45b4700a-399a-42c2-ac5a-07fdbe80195d
# ╠═f8476e50-cdc7-4239-ace4-cfc90fd94e4e
# ╠═033e3abe-318b-4ccc-a18d-f4b692312d92
# ╠═d627fd5a-cddd-4179-8796-4bc967585ead
# ╠═42de57f4-9742-4e13-b28c-ed9b6aa6d8be
# ╠═b13f1b95-68ee-4b24-8584-e7a837ab0b5d
# ╠═0c26b8ce-6ac8-401e-951f-a1fe9f6d272a
