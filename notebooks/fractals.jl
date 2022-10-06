### A Pluto.jl notebook ###
# v0.19.12

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

# ╔═╡ 68fbb48c-22e9-11ed-2e3b-3f933ac78fe0
import Pkg

# ╔═╡ 2ce20e7a-ed7c-4005-96a7-7ab71954037d
begin
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 7cd7c277-86de-49a6-9349-9190be1d24bc
using SparseArrays, LinearAlgebra, Plots, KrylovKit, SpinSymmetry, ElasticArrays, PlutoUI, LsqFit, LaTeXStrings

# ╔═╡ 1879ca84-b3ae-4891-918f-52a593a169d6
TableOfContents()

# ╔═╡ 3116bdd8-0a1f-4682-8524-7af69ffd80eb
html"""<style>main {max-width: 60%;}</style>"""

# ╔═╡ 6c879b09-dc10-4ffd-8fb3-e9f439126c39
begin
	function generate_sierpinski(pos_x::Float64,pos_y::Float64,depth::Int,matrix=false) #unit steps between points
		N = convert(Int,3*(3^depth+1)/2)
		positions = ElasticArray{Float64}(undef,2,3)
		positions[1,1] = pos_x
		positions[2,1] = pos_y
		positions[1,2] = pos_x-0.5*2^depth
		positions[2,2] = pos_y - sqrt(3)/2*2^depth
		positions[1,3] = pos_x+0.5*2^depth
		positions[2,3] = pos_y - sqrt(3)/2*2^depth
		inds = [1,2,3]
		n_set = 3
		if matrix == true
			matrix = zeros(N,N)
			matrix[1,2] = 1
			matrix[2,3] = 1
			matrix[3,1] = 1
			generate_sierpinski!(positions[:,1],positions[:,2],positions[:,3],1,2,3,positions,matrix,depth)
			matrix = matrix + transpose(matrix)
			return positions, matrix
		end
		generate_sierpinski!(positions[:,1],positions[:,2],positions[:,3],positions,depth)
		return positions
	end
	function generate_sierpinski!(A,B,C,positions,depth::Int)
		if depth == 0
			return
		end
		MAB = (A+B)/2
		MBC = (B+C)/2
		MCA = (C+A)/2
		append!(positions,MAB)
		append!(positions,MBC)
		append!(positions,MCA)
		generate_sierpinski!(A,MAB,MCA,positions,depth-1)
		generate_sierpinski!(MAB,B,MBC,positions,depth-1)
		generate_sierpinski!(MCA,MBC,C,positions,depth-1)
	end

	function generate_sierpinski!(A,B,C,indA,indB,indC,positions,matrix,depth::Int)
		if depth == 0
			return
		end
		MAB = (A+B)/2
		MBC = (B+C)/2
		MCA = (C+A)/2
		append!(positions,MAB)
		indAB = size(positions)[2]
		indBC = indAB + 1
		indCA = indBC + 1
		append!(positions,MBC)
		append!(positions,MCA)
		#First make too many connections, reduce later
		matrix[indA,indAB] = 1
		matrix[indB,indAB] = 1
		matrix[indB,indBC] = 1
		matrix[indC,indBC] = 1
		matrix[indC,indCA] = 1
		matrix[indA,indCA] = 1
		matrix[indAB,indBC] = 1
		matrix[indBC,indCA] = 1
		matrix[indCA,indAB] = 1
		generate_sierpinski!(A,MAB,MCA,indA,indAB,indCA,positions,matrix,depth-1)
		generate_sierpinski!(MAB,B,MBC,indAB,indB,indBC,positions,matrix,depth-1)
		generate_sierpinski!(MCA,MBC,C,indCA,indBC,indC,positions,matrix,depth-1)
	end
	
	function reduce_adjacency_matrix!(adj,points) #by construction minimum distance 1
		N = size(points)[1]
		for i in 1:N
			for j in 1:N
				if adj[i,j] == 1
					if !isapprox(norm(points[i,:]-points[j,:]),1.0)
						adj[i,j] = 0
					end
				end
			end
		end
	end
	
end

# ╔═╡ 4fc28aa5-207e-4c48-a78f-71476a4b40de
begin
	function σzσz_1ex(i,j,N)
		res = spdiagm(ones(N))
		res[i,i] = -1
		res[j,j] = -1
		return res
	end
	function construct_H_from_J(J,N::Int,Δ::Float64=0)
		res = spzeros(Float64,N,N)
		for i in 1:N
			for j in i+1:N
				res += J[i,j]*(2*sparse([i],[j],[1],N,N)+2*sparse([j],[i],[1],N,N)+Δ*σzσz_1ex(i,j,N))
			end
		end
		return res
	end
	construct_H_from_J(J) = construct_H_from_J(J,size(J)[1],0.0)
end

# ╔═╡ 440fcfc8-d164-4bde-b45f-95a21d19ca6d
function plot_wavefunction_amplitude(points,matrix,ψ,factor=10,title="")
	probs = abs2.(ψ)
	p = scatter(points[:,1],points[:,2],legend=false,ms=10,alpha=probs*factor,c="red",title=title)
	for i in 1:size(points)[1]		
		for j in 1:i
			if matrix[i,j] == 1
				plot!(p,[points[i,1],points[j,1]],[points[i,2],points[j,2]],c="black",lw=1)
			end
		end
	end
	return p
end

# ╔═╡ e627b35e-86ec-4413-b9a5-098e6b7a771f
function propagate(H,t,ψ)
	return exponentiate(H,im*t,ψ;ishermitian=true)[1]
end

# ╔═╡ 1de1e6ac-0157-41dc-905a-a974d32a3b23
function msd(p,ψ,s=1)
	prob = abs2.(ψ)
	N = length(prob)
	var = 0
	for i in 1:N
		var += sum(abs2.(p[s,:]-p[i,:]))*prob[i]
	end
	return var
end

# ╔═╡ da93d4ec-cfb9-43fb-81d8-585b9aaf8f14
md"# Sierpinski gasket - Wrong connectivity"

# ╔═╡ 4d9a3daa-3da8-4c00-813c-8dcfc8a6bdeb
begin
	depth = 5
	N = convert(Int,3*(3^depth+1)/2)
	p,m = generate_sierpinski(0.,0.,depth,true)
	p = transpose(p)
	reduce_adjacency_matrix!(m,p)
	N
end

# ╔═╡ 37942b14-0087-4c91-92ce-f9af9462224d
maximum(sum(m,dims=1)[:])

# ╔═╡ 01f627d4-9ba6-4402-bb38-074cad101959
begin
	pl = scatter(p[:,1],p[:,2],legend=false,ms=2,axis=nothing)#,texts=1:size(p)[1])
	for i in 1:size(p)[1]		
		for j in 1:i
			if m[i,j] == 1
				plot!(pl,[p[i,1],p[j,1]],[p[i,2],p[j,2]],c="black")
			end
		end
	end
	pl
end

# ╔═╡ 7cd5f604-f06f-46e7-88be-0b479850272e
md"## Spin Transport"

# ╔═╡ 3e52dbae-73dc-442e-8f80-3d0bf9f22353
H = construct_H_from_J(m,size(m)[1],0.)

# ╔═╡ 70e2bf1a-824c-48ca-a392-e5ed05933592
begin
	ψ0 = zeros(N)
	s = 1
	ψ0[s] = 1.
	trange = 0:0.05:1
	NT = length(trange)
	ψ0
end

# ╔═╡ e8046f65-dce6-4bc2-adf3-5da8b2d66544
begin
	ψs = zeros(ComplexF64,length(trange),N)
	msds = zeros(Float64,length(trange))
	for i in 1:length(trange)
		ψs[i,:] = propagate(H,trange[i],ψ0)
		msds[i] = msd(p,ψs[i,:],s)
	end
end

# ╔═╡ 04e56cdf-3314-4b06-9892-ef970eb0df14
ψs

# ╔═╡ 1700ee4b-4608-48d3-9a8f-1d6110ae7410
@bind ti Slider(1:length(trange))

# ╔═╡ 5400efec-3712-4538-b1e2-096b879450e2
plot_wavefunction_amplitude(p,m,ψs[ti,:],100)

# ╔═╡ 5fb6dd50-130d-4f28-93db-22d6aa07bca8
begin
	plot(trange[2:NT],msds[2:NT],xaxis=:log,yaxis=:log,label="MSD")
	plot!(trange[2:NT],trange[2:NT].^2.5 * msds[2]/trange[2]^2.5,label="t^2.5")
	plot!(trange[2:NT],trange[2:NT].^2 * msds[2]/trange[2]^2,label="t^2")
	#plot!(trange[2:NT],trange[2:NT].^1 * msds[2]/trange[2]^1,label="t^1")
	#plot!(trange[2:NT],trange[2:NT].^1.58 * msds[2]/trange[2]^1.58,label="t^1.58")
end

# ╔═╡ 5d1575f6-7303-4b91-a3dc-0f6dfafb46d2
md"## Photon Propagation"

# ╔═╡ 4492a08b-95f8-4447-92f9-5b53a8854c14
begin
	function photonic_H(adj,N::Int,β::Float64,C::Float64)
		res = spzeros(Float64,N,N)
		for i in 1:N
			res[i,i] = β
			for j in 1:N
				if i!=j
					res[i,j] = adj[i,j]*C
				end
			end
		end
		return res
	end
	photonic_H(adj) = photonic_H(adj,size(adj)[1],1.0,1.0)
end

# ╔═╡ f4b22aaa-5b26-45f8-949f-f911a54ef867
H_ph = photonic_H(m,N,0.,5.0)

# ╔═╡ 65120f0c-fa6b-41d3-8d65-9cd8fd08a8c5
H_ph == H

# ╔═╡ c6bf195c-36f4-4803-a813-30d90025c357
ψ0

# ╔═╡ 5f8ecc80-200f-4e73-b101-d20cd8d87fab
begin
	ψs_ph = zeros(ComplexF64,length(trange),N)
	msds_ph = zeros(Float64,length(trange))
	for i in 1:length(trange)
		ψs_ph[i,:] = propagate(H_ph,trange[i],ψ0)
		msds_ph[i] = msd(p,ψs_ph[i,:],s)
	end
end

# ╔═╡ da557dd1-b0e6-4287-a8d3-d76555a4c1e1
@bind ti_ph Slider(1:length(trange))

# ╔═╡ e43d14f9-1154-4b18-a886-168bda51da21
plot_wavefunction_amplitude(p,m,ψs_ph[ti_ph,:],10,trange[ti_ph])

# ╔═╡ 14211d72-ac07-4977-b466-38aa24b52b89
begin
	plot(trange[2:NT],msds_ph[2:NT],xaxis=:log,yaxis=:log,label="MSD")
	plot!(trange[2:NT],trange[2:NT].^2.5 * msds_ph[2]/trange[2]^2.5,label="t^2.5")
	plot!(trange[2:NT],trange[2:NT].^2 * msds_ph[2]/trange[2]^2,label="t^2")
	#plot!(trange[2:NT],trange[2:NT].^1 * msds[2]/trange[2]^1,label="t^1")
	#plot!(trange[2:NT],trange[2:NT].^1.58 * msds[2]/trange[2]^1.58,label="t^1.58")
end

# ╔═╡ 136a1c7c-5c1d-47d1-b8b6-15cd18c97029
md"# Sierpinski gasket - Correct connectivity"

# ╔═╡ 5e888e63-6063-4e88-81bd-d4483ac00868
function nn_connectivity(points) #by construction minimum distance 1
		N = size(points)[1]
		adj = zeros(N,N)
		for i in 1:N
			for j in 1:N
				if isapprox(norm(points[i,:]-points[j,:]),1.0)
					adj[i,j] = 1
				end
			end
		end
		return adj
end

# ╔═╡ 0a80d332-5374-4b4d-8237-da0fc35bfe58
adj = nn_connectivity(p)

# ╔═╡ 274ae4d3-57c8-469c-be55-411a2e18cce8
begin
	plc = scatter(p[:,1],p[:,2],legend=false,ms=2,axis=nothing)#,texts=1:size(p)[1])
	for i in 1:size(p)[1]		
		for j in 1:i
			if adj[i,j] == 1
				plot!(plc,[p[i,1],p[j,1]],[p[i,2],p[j,2]],c="black")
			end
		end
	end
	plc
end

# ╔═╡ 3a6fbbd0-70df-40ce-b81e-3b80b1fa252c
md"## Photon Propagation"

# ╔═╡ 756845f8-6001-4a74-a6b1-383c52d77424
begin
	C = 1.0
	δt = 0.1
	trange_phc = 0:δt:10
	NT_phc = length(trange_phc)
	trange_phc
end

# ╔═╡ dd36f02c-aa49-41a0-9c3d-0c4362d5ff28
H_phc = photonic_H(adj,N,0.,C)

# ╔═╡ e06b1b61-26ae-4a48-aa81-d07e2ed206bc
begin
	ψs_phc = zeros(ComplexF64,length(trange_phc),N)
	msds_phc = zeros(Float64,length(trange_phc))
	for i in 1:length(trange_phc)
		ψs_phc[i,:] = propagate(H_phc,trange_phc[i],ψ0)
		msds_phc[i] = msd(p,ψs_phc[i,:],s)
	end
end

# ╔═╡ 7d8aaf37-609e-4a91-acc0-0b57d9315480
@bind ti_phc Slider(1:length(trange_phc))

# ╔═╡ 6380612f-1ce8-4972-b0a5-a206a6cf6a16
plot_wavefunction_amplitude(p,m,ψs_phc[ti_phc,:],50,trange_phc[ti_phc])

# ╔═╡ 065c6d6f-b659-4e64-8ceb-8747dd9a2440
begin
	plot(trange_phc[2:NT_phc],msds_phc[2:NT_phc],xaxis=:log,yaxis=:log,label="MSD",legend=:bottomright)
	plot!(trange_phc[2:NT_phc],trange_phc[2:NT_phc].^2.5 * msds_phc[2]/trange_phc[2]^2.5,label="t^2.5")
	plot!(trange_phc[2:NT_phc],trange_phc[2:NT_phc].^2 * msds_phc[2]/trange_phc[2]^2,label="t^2")
	plot!(trange_phc[2:NT_phc],trange_phc[2:NT_phc].^1 * msds_phc[2]/trange_phc[2]^1,label="t^1")
	plot!([3*1e-1,3*1e-1],[1e-1,1e1],linestyle=:dash,c="black")
	#plot!(trange[2:NT],trange[2:NT].^1.58 * msds[2]/trange[2]^1.58,label="t^1.58")
end

# ╔═╡ 4ff8d617-22e9-4cab-817c-de288011fc43
begin
	Tcut = 1.2
	k = round(Int,Tcut/δt)
	tdata = 0:δt:Tcut
	rdata = msds_phc[1:k+1]
	@. model(t,params) = params[1] * t^params[2]
	params0 = [1.,1.]
end

# ╔═╡ 7ba173c6-f59b-4d7d-989d-643e0d92be4f
begin
	Tstart = 3
	Tend = 74
	kstart = round(Int,Tstart/δt)
	kend = round(Int,Tend/δt)
	tdata2 = Tstart:δt:Tend
	rdata2 = msds_phc[kstart+1:kend+1]
end

# ╔═╡ ba61f42c-ce9b-4e08-adce-68c0175c9daf
params = curve_fit(model, tdata, rdata, params0).param

# ╔═╡ 5a336720-cf90-4dd6-8fd4-60ff74e84fe3
params2 = curve_fit(model, tdata2, rdata2, params0).param

# ╔═╡ 5f309a3a-68f5-4eaa-954b-ccba296fc37a
begin
	plot(trange_phc[2:NT_phc],msds_phc[2:NT_phc],xaxis=:log,yaxis=:log,label="MSD",legend=nothing,xlabel="Jt [a.u.]",ylabel="MSD [a.u.]")#,xlims=[2*1e-1,3*1e1],ylims=[1e-1,1e3])
	plot!(trange_phc[2:k+1],model(trange_phc[2:k+1],params),ls=:dash,c="black")
	plot!([Tcut,Tcut],[5*1e-2,10*1e3],linestyle=:dashdot,c="black")
	plot!([Tstart,Tstart],[5*1e-2,10*1e3],linestyle=:dashdot,c="black")
	plot!(trange_phc[kstart+1:kend+1],model(trange_phc[kstart+1:kend+1],params2),linestyle=:dash,c="black")
	plot!([Tend,Tend],[5*1e-2,10*1e3],linestyle=:dashdot,c="black")
	annotate!(0.4,0.05,L"\sim(Jt)^{2.5}")
	annotate!(15,100,L"\sim(Jt)^{1.5}")
end

# ╔═╡ 9a94018b-09fb-42e1-9e35-919d4f5d8ec9
md"## Spin Transport"

# ╔═╡ 0665e7b7-4458-43d4-9305-f46ff4bf3484
H_c = construct_H_from_J(adj,size(adj)[1],0.)

# ╔═╡ 3653c03c-7da2-40d3-8fdb-80d33a4150d9
ψ0

# ╔═╡ a31b9a72-04c0-4540-8fc9-f2e0b4ac9760
begin
	ψs_c = zeros(ComplexF64,length(trange),N)
	msds_c = zeros(Float64,length(trange))
	for i in 1:length(trange)
		ψs_c[i,:] = propagate(H_c,trange[i],ψ0)
		msds_c[i] = msd(p,ψs_c[i,:],s)
	end
end

# ╔═╡ 79bd0473-9eac-4b38-82d3-080996a926b7
@bind ti_c Slider(1:length(trange))

# ╔═╡ e31a8d02-3cfb-4e7c-b474-2bcf87615018
plot_wavefunction_amplitude(p,adj,ψs_c[ti_c,:],100,trange[ti_c])

# ╔═╡ 10ce3648-a6bd-4ed4-8c28-85da7fe6d078
begin
	plot(trange[2:NT],msds_c[2:NT],xaxis=:log,yaxis=:log,label="MSD")
	plot!(trange[2:NT],trange[2:NT].^2.5 * msds_c[2]/trange[2]^2.5,label="t^2.5")
	plot!(trange[2:NT],trange[2:NT].^2 * msds_c[2]/trange[2]^2,label="t^2")
	#plot!(trange[2:NT],trange[2:NT].^1 * msds[2]/trange[2]^1,label="t^1")
	#plot!(trange[2:NT],trange[2:NT].^1.58 * msds[2]/trange[2]^1.58,label="t^1.58")
end

# ╔═╡ c3d5a01e-5bca-4ab8-b2e4-cdf64f91b09a


# ╔═╡ Cell order:
# ╠═68fbb48c-22e9-11ed-2e3b-3f933ac78fe0
# ╠═2ce20e7a-ed7c-4005-96a7-7ab71954037d
# ╠═7cd7c277-86de-49a6-9349-9190be1d24bc
# ╠═1879ca84-b3ae-4891-918f-52a593a169d6
# ╠═3116bdd8-0a1f-4682-8524-7af69ffd80eb
# ╠═6c879b09-dc10-4ffd-8fb3-e9f439126c39
# ╠═4fc28aa5-207e-4c48-a78f-71476a4b40de
# ╠═440fcfc8-d164-4bde-b45f-95a21d19ca6d
# ╠═e627b35e-86ec-4413-b9a5-098e6b7a771f
# ╠═1de1e6ac-0157-41dc-905a-a974d32a3b23
# ╠═da93d4ec-cfb9-43fb-81d8-585b9aaf8f14
# ╠═4d9a3daa-3da8-4c00-813c-8dcfc8a6bdeb
# ╠═37942b14-0087-4c91-92ce-f9af9462224d
# ╠═01f627d4-9ba6-4402-bb38-074cad101959
# ╠═7cd5f604-f06f-46e7-88be-0b479850272e
# ╠═3e52dbae-73dc-442e-8f80-3d0bf9f22353
# ╠═70e2bf1a-824c-48ca-a392-e5ed05933592
# ╠═e8046f65-dce6-4bc2-adf3-5da8b2d66544
# ╠═04e56cdf-3314-4b06-9892-ef970eb0df14
# ╠═1700ee4b-4608-48d3-9a8f-1d6110ae7410
# ╠═5400efec-3712-4538-b1e2-096b879450e2
# ╠═5fb6dd50-130d-4f28-93db-22d6aa07bca8
# ╠═5d1575f6-7303-4b91-a3dc-0f6dfafb46d2
# ╠═4492a08b-95f8-4447-92f9-5b53a8854c14
# ╠═f4b22aaa-5b26-45f8-949f-f911a54ef867
# ╠═65120f0c-fa6b-41d3-8d65-9cd8fd08a8c5
# ╠═c6bf195c-36f4-4803-a813-30d90025c357
# ╠═5f8ecc80-200f-4e73-b101-d20cd8d87fab
# ╠═da557dd1-b0e6-4287-a8d3-d76555a4c1e1
# ╠═e43d14f9-1154-4b18-a886-168bda51da21
# ╠═14211d72-ac07-4977-b466-38aa24b52b89
# ╠═136a1c7c-5c1d-47d1-b8b6-15cd18c97029
# ╠═5e888e63-6063-4e88-81bd-d4483ac00868
# ╠═0a80d332-5374-4b4d-8237-da0fc35bfe58
# ╠═274ae4d3-57c8-469c-be55-411a2e18cce8
# ╠═3a6fbbd0-70df-40ce-b81e-3b80b1fa252c
# ╠═756845f8-6001-4a74-a6b1-383c52d77424
# ╠═dd36f02c-aa49-41a0-9c3d-0c4362d5ff28
# ╠═e06b1b61-26ae-4a48-aa81-d07e2ed206bc
# ╠═7d8aaf37-609e-4a91-acc0-0b57d9315480
# ╠═6380612f-1ce8-4972-b0a5-a206a6cf6a16
# ╠═065c6d6f-b659-4e64-8ceb-8747dd9a2440
# ╠═4ff8d617-22e9-4cab-817c-de288011fc43
# ╠═7ba173c6-f59b-4d7d-989d-643e0d92be4f
# ╠═ba61f42c-ce9b-4e08-adce-68c0175c9daf
# ╠═5a336720-cf90-4dd6-8fd4-60ff74e84fe3
# ╠═5f309a3a-68f5-4eaa-954b-ccba296fc37a
# ╠═9a94018b-09fb-42e1-9e35-919d4f5d8ec9
# ╠═0665e7b7-4458-43d4-9305-f46ff4bf3484
# ╠═3653c03c-7da2-40d3-8fdb-80d33a4150d9
# ╠═a31b9a72-04c0-4540-8fc9-f2e0b4ac9760
# ╠═79bd0473-9eac-4b38-82d3-080996a926b7
# ╠═e31a8d02-3cfb-4e7c-b474-2bcf87615018
# ╠═10ce3648-a6bd-4ed4-8c28-85da7fe6d078
# ╠═c3d5a01e-5bca-4ab8-b2e4-cdf64f91b09a
