### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ 02d2b47f-044e-4715-82f7-5d7a062d17a2
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 14a93042-910e-11ed-03e3-7b58181ed649
using LinearAlgebra,Plots,JLD2,Statistics,PlutoUI

# ╔═╡ 4577466d-7fdf-4cbf-b300-decdc7d2ee06
TableOfContents()

# ╔═╡ dd1cd320-f66a-4a66-9814-274aad25d324
T = 0.1:0.1:20.0

# ╔═╡ 90383e3f-9d72-4447-8c65-a7738ee90b79
path = pwd() * "/data/with_otoc"

# ╔═╡ e48aa741-b660-4145-9af4-ebf02f92eb87
datapaths = [path*"/WAHUHA/N = 11/",path*"/WAHUHA/N = 12/",path*"/WAHUHA/N = 13/",
path*"/Rhim71/N = 11/",path*"/Rhim71/N = 12/",path*"/Rhim71/N = 13/",
path*"/WAHUHA_FR/N = 11/",path*"/WAHUHA_FR/N = 12/",path*"/WAHUHA_FR/N = 13/",
path*"/Rhim71_FR/N = 11/",path*"/Rhim71_FR/N = 12/",path*"/Rhim71_FR/N = 13/"
]

# ╔═╡ 7a514d8c-38d3-4d60-8efb-ed8463b615a1
begin
	function state_mean(A)
		return mean(A,;dims=4)[:,:,:,1]
	end
	function reduce(A)
		return A[:,:,1]
	end
end

# ╔═╡ 31fcd6a6-0a91-4610-8f0f-e0327f1873d3
md"# N = 11"

# ╔═╡ ab930454-8faf-4f22-a923-fe6c0dce68d8
begin
	fids1 = load(datapaths[1]*"468288_N11_BS.jld2","fidelities")
	fids5 = load(datapaths[1]*"468289_N11_BS.jld2","fidelities")
	fids10 = load(datapaths[1]*"468290_N11_BS.jld2","fidelities")
	fids50 = load(datapaths[1]*"468291_N11_BS.jld2","fidelities")
	fids100 = load(datapaths[1]*"468292_N11_BS.jld2","fidelities")
	fids500 = load(datapaths[1]*"468293_N11_BS.jld2","fidelities")
	fids1000 = load(datapaths[1]*"468294_N11_BS.jld2","fidelities")
	
	otocs1 = load(datapaths[1]*"468288_N11_BS.jld2","otocs")
	otocs5 = load(datapaths[1]*"468289_N11_BS.jld2","otocs")
	otocs10 = load(datapaths[1]*"468290_N11_BS.jld2","otocs")
	otocs50 = load(datapaths[1]*"468291_N11_BS.jld2","otocs")
	otocs100 = load(datapaths[1]*"468292_N11_BS.jld2","otocs")
	otocs500 = load(datapaths[1]*"468293_N11_BS.jld2","otocs")
	otocs1000 = load(datapaths[1]*"468294_N11_BS.jld2","otocs")
end

# ╔═╡ cd44e8fd-3efc-467b-b2ec-9ef0c1623dcb
reduce(state_mean(otocs10))

# ╔═╡ 10b781e2-10bf-4c94-83a2-377b7737499f
size(fids1)

# ╔═╡ 3700073a-b14b-455f-aa1b-68e828d7ad18
begin
	plot(T,fids500[:,1,1],legend=nothing)
	for i in 2:49
		plot!(T,fids500[:,1,i])
	end
	plot!(T,fids500[:,1,50])
end

# ╔═╡ 82efb7b0-9595-4062-82d3-e8811b23ff96


# ╔═╡ 11877b7b-8ee3-430e-8f4e-5d26931d735a
size(otocs100)

# ╔═╡ ac1b47ae-567c-4440-85c8-9e688f0447d0
begin
	plot(T,otocs100[:,3,1,1],legend=nothing)
	plot!(T,reduce(state_mean(otocs100))[:,3])
end

# ╔═╡ d21444a3-188d-445b-b11a-b55fb961109f
begin
	plot(T,otocs100[:,1,1,1],legend=nothing)
	for i in 2:49
		plot!(T,otocs100[:,1,1,i])
	end
	plot!(T,otocs100[:,1,1,50])
end

# ╔═╡ Cell order:
# ╠═14a93042-910e-11ed-03e3-7b58181ed649
# ╠═02d2b47f-044e-4715-82f7-5d7a062d17a2
# ╠═4577466d-7fdf-4cbf-b300-decdc7d2ee06
# ╠═dd1cd320-f66a-4a66-9814-274aad25d324
# ╠═90383e3f-9d72-4447-8c65-a7738ee90b79
# ╠═e48aa741-b660-4145-9af4-ebf02f92eb87
# ╠═7a514d8c-38d3-4d60-8efb-ed8463b615a1
# ╠═31fcd6a6-0a91-4610-8f0f-e0327f1873d3
# ╠═ab930454-8faf-4f22-a923-fe6c0dce68d8
# ╠═cd44e8fd-3efc-467b-b2ec-9ef0c1623dcb
# ╠═10b781e2-10bf-4c94-83a2-377b7737499f
# ╠═3700073a-b14b-455f-aa1b-68e828d7ad18
# ╠═82efb7b0-9595-4062-82d3-e8811b23ff96
# ╠═11877b7b-8ee3-430e-8f4e-5d26931d735a
# ╠═ac1b47ae-567c-4440-85c8-9e688f0447d0
# ╠═d21444a3-188d-445b-b11a-b55fb961109f
