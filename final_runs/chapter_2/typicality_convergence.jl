### A Pluto.jl notebook ###
# v0.19.17

using Markdown
using InteractiveUtils

# ╔═╡ 079ca122-8c2c-4e3d-9eee-d936361371b5
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 1fa345ce-4a5e-11ed-2fe2-75fc302d3420
using LinearAlgebra,JLD2,Statistics,PlutoUI, SpinSymmetry, BenchmarkTools

# ╔═╡ 54a9508a-7182-42cd-870e-782d1f57ebb1
using Plots

# ╔═╡ c37982e1-313e-4223-a5a8-580a689a1e3a
TableOfContents()

# ╔═╡ 383c87ea-5715-4f1a-906c-9874d2e66ada
begin
	function state_mean(A,n_states)
			return mean(A[:,:,:,1:n_states];dims=4)[:,:,:,1]
		end
		
		function state_std(A,n_states)
			return std(A[:,:,:,1:n_states];dims=4)[:,:,:,1]
		end
end

# ╔═╡ 2fb825be-cb30-477e-b11c-d13c895c1af4
begin
	δt = 0.1
	T = 10.0
	Ns = [5,6,7,8,9,10,11,12,13,14,15,16,17]
	trange = 0.0:δt:T
	ts = length(trange)
end

# ╔═╡ 60f87fc0-f3a0-4a43-96b1-0532fc00b672
begin
	dims = Vector{Int64}(undef,13)
	for (i,n) in enumerate(Ns)
		k = div(n-1,2)+1
		dims[i] = basissize(symmetrized_basis(n,k))
	end
end

# ╔═╡ 29501d68-aa87-4a2e-82a5-4287d87459d3
dims

# ╔═╡ 42bf99ad-a834-47a8-9400-9eb1a17d2fbe
md"# Load data"

# ╔═╡ cd98b70a-9061-4502-98c8-e70fcd44baa9
begin
	params1 = Vector{SimulationParamsED}(undef,13)
	params2 = Vector{SimulationParamsED}(undef,13)
	params3 = Vector{SimulationParams}(undef,13)
	params4 = Vector{SimulationParamsED}(undef,13)
	params5 = Vector{SimulationParamsED}(undef,11)
	params6 = Vector{SimulationParams}(undef,11)
	
	
	params = [params1,params2,params3,params4,params5,params6]
end

# ╔═╡ a69b0740-8fd7-4cfb-bad5-c90d5ef163a5
begin
	data_EDtr_sector = Vector{Array{Float64,3}}(undef,13)
	data_ED_sector = Vector{Array{Float64,4}}(undef,13)
	data_Krylov_sector = Vector{Array{Float64,4}}(undef,13)
	data_EDtr_total = Vector{Array{Float64,3}}(undef,11)
	data_ED_total = Vector{Array{Float64,4}}(undef,11)
	data_Krylov_total = Vector{Array{Float64,4}}(undef,11)
	data = [data_EDtr_sector,data_ED_sector,data_Krylov_sector,data_EDtr_total,data_ED_total,data_Krylov_total]
end

# ╔═╡ c8459c23-67c0-4abe-81c3-5033b8374450
path = pwd()*"/data/"

# ╔═╡ aa28d126-83b4-4af9-8a30-1763eabaf3d3
folders = ["ED_tr/sector/","ED/sector/","Krylov/sector/","ED_tr/total/","ED/total/","Krylov/total/"]

# ╔═╡ 8d4a6591-7425-4dcb-80e4-55fa60f6cb42
load(path*folders[1]*"7333905_N6_ED.jld2","params")

# ╔═╡ 1670e4ea-4488-455f-840d-5d95038c3cb2
begin
	tmp = 1
	tmpdata = Vector{Array{Float64,3}}(undef,3)
	for (i,folder) in enumerate(folders)
			print(folder,"\n")
			for (j,filename) in enumerate(readdir(path*folders[i]))
				if !occursin("otoc",filename)
					print(j,":",filename,"\n")
					if occursin("_tr",folder)
						if occursin("sector",folder)
							factor = dims[j]
						else
							factor = 2^Ns[j]
						end
						if occursin("sector",folder)
							if occursin("N16",filename) || occursin("N17",filename)
								tlen = 11
							else
								tlen = 101
							end
							data[i][j] = 2*ones(tlen,Ns[j],1)-2*load(path*folder*filename,"data")./factor
						else
							if occursin("N15",filename)
								if tmp < 3
									tlen = 4
								else
									tlen = 3
								end
								
								tmpdata[tmp] = 2*ones(tlen,15,1)-2*load(path*folder*filename,"data")./factor
								if tmp == 3
									print(tmp)
									data[i][11] = vcat(tmpdata[1],tmpdata[2],tmpdata[3])
								end
								tmp += 1
							elseif occursin("N14",filename)
								tlen = 11
								data[i][j] = 2*ones(tlen,Ns[j],1)-2*load(path*folder*filename,"data")./factor
							else
								tlen = 101
								data[i][j] = 2*ones(tlen,Ns[j],1)-2*load(path*folder*filename,"data")./factor
							end
						end		
					else
						if occursin("sector",folder) && occursin("ED",folder)
							if occursin("N16",filename)
								tlen = 31
								nstates = 10
							elseif occursin("N17",filename)
								tlen = 11
								nstates = 10
							else
								tlen = 101
								nstates = 1000
							end
						elseif occursin("sector",folder) && occursin("Krylov",folder)
							tlen = 101
							if occursin("N16",filename) || occursin("N17",filename)
								nstates = 10
							else
								nstates = 1000
							end
						elseif occursin("total",folder) && occursin("ED",folder)
							if occursin("N14",filename)
								tlen = 31
								nstates = 10
							elseif occursin("N15",filename)
								tlen = 11
								nstates = 10
							else
								tlen = 101
								nstates = 1000
							end
						elseif occursin("total",folder) && occursin("Krylov",folder)
							if occursin("N14",filename) || occursin("N15",filename)
								tlen = 51
								nstates = 10
							else
								tlen = 101
								nstates = 1000
							end
						end
						data[i][j] = 2*ones(tlen,Ns[j],1,nstates)-2*load(path*folder*filename,"data")
					end
					params[i][j] = load(path*folder*filename,"params")
				end
			end
		end
end

# ╔═╡ 51158cf8-e121-499a-b36e-29fb8e663155
md"# System Size vs. Error"

# ╔═╡ 79a71641-0fe2-4777-8673-e1cc8b539be1
Ns[5]

# ╔═╡ 04d4af1f-f2ef-47ba-b5e1-33cb93c3faf8
begin
	k = 6
	plot(0:0.1:1,state_mean(data_ED_total[k],10)[1:11,:,1])
	plot!(0:0.1:1,data_EDtr_total[k][1:11,:,1],legend=nothing,color="black")
end

# ╔═╡ 3e3a86aa-1600-4337-b448-480c59bbcfc9
begin
	k2 = 11
	plot(0:0.1:1,state_mean(data_ED_total[k2],10)[1:11,:,1,1])
	plot!(0:0.1:1,data_EDtr_total[k2][1:11,:,1],legend=nothing,color="black")
end

# ╔═╡ 031cb4f0-8628-49ca-bbf2-73f9f53dbf3b
(2*ones(ts,5,1)-data_EDtr_sector[1])*10

# ╔═╡ cdcc9475-6305-48eb-9376-ebd0b4f1c8e8
plot(trange,data_EDtr_total[9][:,:,1])

# ╔═╡ 523a0091-acfa-41f1-be50-6d8d753fda3a
2*0.01

# ╔═╡ 03abf2d2-3461-4e19-822f-1428304b008c
plot(abs.(data_EDtr_total[9][2:ts,:,1]-state_mean(data_Krylov_total[9],10)[2:ts,:])./data_EDtr_total[9][2:ts,:,1]*100)

# ╔═╡ 6222f35e-b931-41d1-a847-89d9f0ce1cc8
data_EDtr_total[9][2:21,:,1]

# ╔═╡ f7e68416-5ccb-49af-a80f-646e2af64ac2
maximum(abs.(data_EDtr_total[8][2:ts,:,1]-state_mean(data_ED_total[8],10)[2:ts,:])./data_EDtr_total[8][2:ts,:,1]*100)

# ╔═╡ 31b41f61-95bf-4ba4-abdd-927dcd2e8045
begin
	errors_ED_sector = Vector{Float64}(undef,13)
	errors_Krylov_sector = Vector{Float64}(undef,13)
	
	errors_ED_total = Vector{Float64}(undef,11)
	errors_Krylov_total = Vector{Float64}(undef,11)

	stdEs = Vector{Float64}(undef,13)
	stdKs = Vector{Float64}(undef,13)

	stdEt = Vector{Float64}(undef,11)
	stdKt = Vector{Float64}(undef,11)

	states = 10
	Tmax = 21
	for i in 1:11
		errors_ED_sector[i] = mean(abs.(data_EDtr_sector[i][1:Tmax,:,1]-state_mean(data_ED_sector[i],states)[1:Tmax,:,1]))
		
		errors_Krylov_sector[i] = mean(abs.(data_EDtr_sector[i][1:Tmax,:,1]-state_mean(data_Krylov_sector[i],states)[1:Tmax,:,1]))

		stdEs[i] = std(mean(abs.(data_EDtr_sector[i][1:Tmax,:,1]-data_ED_sector[i][1:Tmax,:,1,j])) for j in 1:states)/sqrt(states)

		stdKs[i] = std(mean(abs.(data_EDtr_sector[i][1:Tmax,:,1]-data_Krylov_sector[i][1:Tmax,:,1,j])) for j in 1:states)/sqrt(states)
	end
	for i in 1:9
		errors_ED_total[i] = mean(abs.(data_EDtr_total[i][1:Tmax,:,1]-state_mean(data_ED_total[i],states)[1:Tmax,:,1]))
		
		errors_Krylov_total[i] = mean(abs.(data_EDtr_total[i][1:Tmax,:,1]-state_mean(data_Krylov_total[i],states)[1:Tmax,:,1]))

		stdEt[i] = std(mean(abs.(data_EDtr_total[i][1:Tmax,:,1]-data_ED_total[i][1:Tmax,:,1,j])) for j in 1:states)/sqrt(states)

		stdKt[i] = std(mean(abs.(data_EDtr_total[i][1:Tmax,:,1]-data_Krylov_total[i][1:Tmax,:,1,j])) for j in 1:states)/sqrt(states)
	end
end

# ╔═╡ 73c95dc6-85ec-4c1d-b720-d242c987af4b
begin
	plot(Ns[1:11],errors_ED_sector[1:11],yerror=stdEs[1:11],label="ED sector",xlabel="N",ylabel="ϵ",xticks=Ns,yminorticks=true)#,ylim=[1e-3,2e-1])
	plot!(Ns[1:11],errors_Krylov_sector[1:11],yerror=stdKs[1:11],label="Krylov sector")

	plot!(Ns[1:9],errors_ED_total[1:9],yerror=stdEt[1:9],label="ED total")
	plot!(Ns[1:9],errors_Krylov_total[1:9],yerror=stdKt[1:9],label="Krylov total",yaxis=:log)
end

# ╔═╡ ef1518f2-ab35-4c98-bdfd-01397f475594
jldopen("errors.jld2", "w") do file
    file["ED_sec"] = errors_ED_sector
	file["stdED_sec"] = stdEs
	file["Kr_sec"] = errors_Krylov_sector
	file["stdKr_sec"] = stdKs
	file["ED_tot"] = errors_ED_total
	file["stdED_tot"] = stdEt
	file["Kr_tot"] = errors_Krylov_total
	file["stdKr_tot"] = stdKt
    end

# ╔═╡ b91a2315-75d8-4086-b85e-aba105e9e305
jldopen("data.jld2", "w") do file
    file["tr"] = data_EDtr_total[9]
	file["ED"] = data_ED_total[9]
	file["Kr"] = data_Krylov_total[9]
    end

# ╔═╡ 82b08be4-f57d-4cf1-87d9-c1b6b6516466
size(data_ED_total[9])

# ╔═╡ fedfb82f-805d-4992-83fc-df4c4c937fba
md"# System Size vs. Error - Multiple states"

# ╔═╡ b373574b-6d45-4d07-8563-6231a597beea
# ╠═╡ disabled = true
#=╠═╡
begin
	mean_ED_sector = Vector{Array{Float64,3}}(undef,11)
	mean_Krylov_sector = Vector{Array{Float64,3}}(undef,11)
	mean_ED_total = Vector{Array{Float64,3}}(undef,9)
	mean_Krylov_total = Vector{Array{Float64,3}}(undef,9)

	std_ED_sector = Vector{Array{Float64,3}}(undef,11)
	std_Krylov_sector = Vector{Array{Float64,3}}(undef,11)
	std_ED_total = Vector{Array{Float64,3}}(undef,9)
	std_Krylov_total = Vector{Array{Float64,3}}(undef,9)

	errors2_ED_sector = Matrix{Float64}(undef,11,999)
	errors2_Krylov_sector = Matrix{Float64}(undef,11,999)
	errors2_ED_total = Matrix{Float64}(undef,9,999)
	errors2_Krylov_total = Matrix{Float64}(undef,9,999)
end
  ╠═╡ =#

# ╔═╡ c6cf9fc9-f06e-466b-890b-fe3291a1610b
# ╠═╡ disabled = true
#=╠═╡
begin
	for i in 1:11
		mean_ED_sector[i] = zeros(ts,Ns[i],999)
		mean_Krylov_sector[i] = zeros(ts,Ns[i],999)
		std_ED_sector[i] = zeros(ts,Ns[i],999)
		std_Krylov_sector[i] = zeros(ts,Ns[i],999)
		for s in 1:999
		mean_ED_sector[i][:,:,s] = state_mean(data_ED_sector[i],s+1)[:,:,1]
		mean_Krylov_sector[i][:,:,s] = state_mean(data_ED_sector[i],s+1)[:,:,1]
		std_ED_sector[i][:,:,s] = state_std(data_ED_sector[i],s+1)[:,:,1]./sqrt(s)
		std_Krylov_sector[i][:,:,s] = state_std(data_ED_sector[i],s+1)[:,:,1]./sqrt(s)
		
		errors2_ED_sector[i,s] = mean(abs.(data_EDtr_sector[i][:,:,1]-mean_ED_sector[i][:,:,s]))
		errors2_Krylov_sector[i,s] = mean(abs.(data_EDtr_sector[i][:,:,1]-mean_Krylov_sector[i][:,:,s]))
		end
	end
	
	for i in 1:9
		mean_ED_total[i] = zeros(ts,Ns[i],999)
		mean_Krylov_total[i] = zeros(ts,Ns[i],999)
		std_ED_total[i] = zeros(ts,Ns[i],999)
		std_Krylov_total[i] = zeros(ts,Ns[i],999)
		for s in 1:999
		mean_ED_total[i][:,:,s] = state_mean(data_ED_total[i],s+1)[:,:,1]
		mean_Krylov_total[i][:,:,s] = state_mean(data_ED_total[i],s+1)[:,:,1]
		std_ED_total[i][:,:,s] = state_std(data_ED_total[i],s+1)[:,:,1]./sqrt(s)
		std_Krylov_total[i][:,:,s] = state_std(data_ED_total[i],s+1)[:,:,1]./sqrt(s)

		errors2_ED_total[i,s] = mean(abs.(data_EDtr_total[i][:,:,1]-mean_ED_total[i][:,:,s]))
		errors2_Krylov_total[i,s] = mean(abs.(data_EDtr_total[i][:,:,1]-mean_Krylov_total[i][:,:,s]))
		end
	end
end
  ╠═╡ =#

# ╔═╡ a0da8248-9838-4207-b047-fd4ce49f5a59
# ╠═╡ disabled = true
#=╠═╡
begin
	plot(2:1000,errors2_ED_sector[1,:],yaxis=:log,xaxis=:log,label="N = 5",xlabel="number of sample states")
	plot!(2:1000,errors2_ED_sector[2,:],yaxis=:log,xaxis=:log,label="N = 6")
	plot!(2:1000,errors2_ED_sector[3,:],yaxis=:log,xaxis=:log,label="N = 7")
	plot!(2:1000,errors2_ED_sector[4,:],yaxis=:log,xaxis=:log,label="N = 8")
	plot!(2:1000,errors2_ED_sector[5,:],yaxis=:log,xaxis=:log,label="N = 9")
	plot!(2:1000,errors2_ED_sector[6,:],yaxis=:log,xaxis=:log,label="N = 10")
	plot!(2:1000,errors2_ED_sector[7,:],yaxis=:log,xaxis=:log,label="N = 11")
	plot!(2:1000,errors2_ED_sector[8,:],yaxis=:log,xaxis=:log,label="N = 12")
	plot!(2:1000,errors2_ED_sector[9,:],yaxis=:log,xaxis=:log,label="N = 13")
	plot!(2:1000,errors2_ED_sector[10,:],yaxis=:log,xaxis=:log,label="N = 14")
	plot!(2:1000,errors2_ED_sector[11,:],yaxis=:log,xaxis=:log,label="N = 15")
end
  ╠═╡ =#

# ╔═╡ eab78ce6-fdba-4e4e-ac21-646ce1ba434b
# ╠═╡ disabled = true
#=╠═╡
begin
	plot(2:1000,errors2_ED_total[1,:],yaxis=:log,xaxis=:log,label="N = 5",xlabel="number of sample states")
	plot!(2:1000,errors2_ED_total[2,:],yaxis=:log,xaxis=:log,label="N = 6")
	plot!(2:1000,errors2_ED_total[3,:],yaxis=:log,xaxis=:log,label="N = 7")
	plot!(2:1000,errors2_ED_total[4,:],yaxis=:log,xaxis=:log,label="N = 8")
	plot!(2:1000,errors2_ED_total[5,:],yaxis=:log,xaxis=:log,label="N = 9")
	plot!(2:1000,errors2_ED_total[6,:],yaxis=:log,xaxis=:log,label="N = 10")
	plot!(2:1000,errors2_ED_total[7,:],yaxis=:log,xaxis=:log,label="N = 11")
	plot!(2:1000,errors2_ED_total[8,:],yaxis=:log,xaxis=:log,label="N = 12")
	plot!(2:1000,errors2_ED_total[9,:],yaxis=:log,xaxis=:log,label="N = 13")
end
  ╠═╡ =#

# ╔═╡ cfeaf277-185c-45bd-b91f-4a77613c777e


# ╔═╡ 405d95ff-bf45-44f2-ba12-301078f22a84


# ╔═╡ c385d15e-1e2f-4732-9173-6d217d2aaf37
# ╠═╡ disabled = true
#=╠═╡
plot([2^9,2^10,2^11,2^12,2^13],[71,535,4455,31297,221948],yaxis=:log,xaxis=:log)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═1fa345ce-4a5e-11ed-2fe2-75fc302d3420
# ╠═079ca122-8c2c-4e3d-9eee-d936361371b5
# ╠═54a9508a-7182-42cd-870e-782d1f57ebb1
# ╠═c37982e1-313e-4223-a5a8-580a689a1e3a
# ╠═383c87ea-5715-4f1a-906c-9874d2e66ada
# ╠═2fb825be-cb30-477e-b11c-d13c895c1af4
# ╠═60f87fc0-f3a0-4a43-96b1-0532fc00b672
# ╠═29501d68-aa87-4a2e-82a5-4287d87459d3
# ╠═42bf99ad-a834-47a8-9400-9eb1a17d2fbe
# ╠═8d4a6591-7425-4dcb-80e4-55fa60f6cb42
# ╠═cd98b70a-9061-4502-98c8-e70fcd44baa9
# ╠═a69b0740-8fd7-4cfb-bad5-c90d5ef163a5
# ╠═c8459c23-67c0-4abe-81c3-5033b8374450
# ╠═aa28d126-83b4-4af9-8a30-1763eabaf3d3
# ╠═1670e4ea-4488-455f-840d-5d95038c3cb2
# ╠═51158cf8-e121-499a-b36e-29fb8e663155
# ╠═79a71641-0fe2-4777-8673-e1cc8b539be1
# ╠═04d4af1f-f2ef-47ba-b5e1-33cb93c3faf8
# ╠═3e3a86aa-1600-4337-b448-480c59bbcfc9
# ╠═031cb4f0-8628-49ca-bbf2-73f9f53dbf3b
# ╠═cdcc9475-6305-48eb-9376-ebd0b4f1c8e8
# ╠═523a0091-acfa-41f1-be50-6d8d753fda3a
# ╠═03abf2d2-3461-4e19-822f-1428304b008c
# ╠═6222f35e-b931-41d1-a847-89d9f0ce1cc8
# ╠═f7e68416-5ccb-49af-a80f-646e2af64ac2
# ╠═31b41f61-95bf-4ba4-abdd-927dcd2e8045
# ╠═73c95dc6-85ec-4c1d-b720-d242c987af4b
# ╠═ef1518f2-ab35-4c98-bdfd-01397f475594
# ╠═b91a2315-75d8-4086-b85e-aba105e9e305
# ╠═82b08be4-f57d-4cf1-87d9-c1b6b6516466
# ╠═fedfb82f-805d-4992-83fc-df4c4c937fba
# ╠═b373574b-6d45-4d07-8563-6231a597beea
# ╠═c6cf9fc9-f06e-466b-890b-fe3291a1610b
# ╠═a0da8248-9838-4207-b047-fd4ce49f5a59
# ╠═eab78ce6-fdba-4e4e-ac21-646ce1ba434b
# ╠═cfeaf277-185c-45bd-b91f-4a77613c777e
# ╠═405d95ff-bf45-44f2-ba12-301078f22a84
# ╠═c385d15e-1e2f-4732-9173-6d217d2aaf37
