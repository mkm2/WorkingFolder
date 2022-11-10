### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 66e8d6f5-7bd2-4894-a9e7-823c9be8d60b
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ df817bae-9639-44aa-bfdd-ecbc5c6b6f7e
using LinearAlgebra,JLD2,Statistics,PlutoUI, SpinSymmetry, BenchmarkTools

# ╔═╡ 97f67aaa-c3d4-4609-ade1-e92943b02c9c
using Plots

# ╔═╡ 7af4ead4-0c87-497e-b100-566f3695ede4
begin
	function disorder_mean(A,n_shots)
				return mean(A[:,:,1:n_shots];dims=3)[:,:,1]
	end
			
	function disorder_std(A,n_shots)
				return std(A[:,:,1:n_shots];dims=3)[:,:,1]
	end
end

# ╔═╡ 12ce7fb1-ad59-4094-ad74-c3f356e6b3be
path = pwd()*"/"

# ╔═╡ 90bfce21-cefe-456c-8865-9f5e586f2023
trange = 10. .^LinRange(-3,8,100)

# ╔═╡ e00ef8bb-d5e3-4459-a214-c103c26a27a8
files = ["7353032_N13_ED.jld2","7361030_N13_ED.jld2","7361031_N13_ED.jld2","7361032_N13_ED.jld2","7361033_N13_ED.jld2","7361034_N13_ED.jld2","7353038_N13_ED.jld2","7361035_N13_ED.jld2","7353040_N13_ED.jld2","7353041_N13_ED.jld2","7361037_N13_ED.jld2","7361038_N13_ED.jld2"]

# ╔═╡ 67a4e20e-3198-4a61-be76-45c8b11265d0
begin
	disorder = Vector{Float64}(undef,12)
	data = Vector{Array{Float64,3}}(undef,12)
	for (i,f) in enumerate(files)
		disorder[i] = load(path*f,"params").DISORDER_PARAM
		data[i] = max.(1e-13,2*ones(100,13,50)-2*load(path*f,"data"))
	end
end

# ╔═╡ b841e087-edb1-41f0-91f8-5a9bd84367cc
disorder

# ╔═╡ 8d58599f-b79d-4160-9c8d-869554430360
begin
	C6 = 400/(2*π) * 1e9 #Hz(μm)^6
	r = 7 #μm
	J = C6/r^6 *0.25 #0.25 from spins
end

# ╔═╡ 9414a48c-fe72-4c07-80fc-01d1d78c2e73
J

# ╔═╡ 4e203efc-bea7-47d7-ace5-50da4f43115b
trange_J = trange/J * 1e6

# ╔═╡ 174e34f7-b948-4338-a272-5b69d65cfe50
begin
	otocs = Vector{Array{Float64,2}}(undef,12)
	otocs_err = Vector{Array{Float64,2}}(undef,12)
	for i in 1:length(disorder)
		otocs[i] = disorder_mean(data[i],50)
		otocs_err[i] = disorder_std(data[i],50)./sqrt(50)
	end
end

# ╔═╡ 0a8e9c9e-185f-431b-8e8f-89eea9c17ef4
(2*rand()-1)

# ╔═╡ 0e9d1f52-837b-445e-9d2c-cebbb8fcdd05
plot(trange_J*6,otocs[1],ribbon=otocs_err[1],xaxis=:log,title="H = 0", legend=nothing,yaxis=:linear,xlabel="t [μs]",ylabel="OTOC",xlim=[1e-1,1e3])

# ╔═╡ bedfd8a7-59cb-4906-9a4e-220e5a0a1a6c
plot(trange_J,otocs[2],ribbon=otocs_err[2],xaxis=:log,title="H = 1", legend=nothing,yaxis=:linear,xlim=[1e-3,1e1],xlabel="t [μs]",ylabel="OTOC")

# ╔═╡ dba9941e-609c-4852-a921-1820ccb678da
plot(trange_J,otocs[3],ribbon=otocs_err[3],xaxis=:log,title="H = 2", legend=nothing,xlabel="t [μs]",ylabel="OTOC",xlim=[1e-3,1e1])

# ╔═╡ ba99f89f-527c-478e-880d-b99fd71fa67e
plot(trange_J,otocs[4],ribbon=otocs_err[4],xaxis=:log,title="H = 3", legend=nothing,yaxis=:linear,xlabel="t [μs]",ylabel="OTOC",xlim=[1e-3,1e1])

# ╔═╡ 0913e415-9308-4ad5-a16e-cbda2b6721a6
plot(trange_J,otocs[5],ribbon=otocs_err[5],xaxis=:log,title="H = 4", legend=nothing,xlabel="t [μs]",ylabel="OTOC",xlim=[1e-3,1e1])

# ╔═╡ 98acf4a3-7a07-4039-9772-350622a3bb9a
plot(trange_J,otocs[6],ribbon=otocs_err[6],xaxis=:log,title="H = 5", legend=nothing,xlabel="t [μs]",ylabel="OTOC",xlim=[1e-3,1e1])

# ╔═╡ 99cd41f2-1d32-4719-9b69-093623e5fdd4
plot(trange_J,otocs[7],ribbon=otocs_err[7],xaxis=:log,title="H = 6", legend=nothing,xlabel="t [μs]",ylabel="OTOC",xlim=[1e-3,1e1])

# ╔═╡ 2fb06cd2-d0d6-46c7-bb83-544435286c51
plot(trange_J,otocs[8],ribbon=otocs_err[8],xaxis=:log,title="H = 7", legend=nothing,xlabel="t [μs]",ylabel="OTOC",xlim=[1e-3,1e1])

# ╔═╡ 5c6536dd-78b9-4403-82f3-3b340916d805
plot(trange_J,otocs[9],ribbon=otocs_err[9],xaxis=:log,title="H = 8", legend=nothing,xlabel="t [μs]",ylabel="OTOC",xlim=[1e-3,1e1])

# ╔═╡ 4d4ad00e-0a23-4f70-bbc1-de17c278377e
plot(trange_J,otocs[10],ribbon=otocs_err[10],xaxis=:log,title="H = 9", legend=nothing,xlabel="t [μs]",ylabel="OTOC",xlim=[1e-3,1e1])

# ╔═╡ e9554486-65cd-4a09-b2cc-1853ce075c00
plot(trange_J,otocs[11],ribbon=otocs_err[11],xaxis=:log,title="H = 11", legend=nothing,xlabel="t [μs]",ylabel="OTOC",xlim=[1e-3,1e1])

# ╔═╡ 7173829a-ca92-4639-bc9c-cf545ebaa42a
plot(trange_J*6,otocs[12],ribbon=otocs_err[12],xaxis=:log,title="H = 12", legend=nothing,xlabel="t [μs]",ylabel="OTOC",xlim=[1e-1,1e3])

# ╔═╡ e9ca6017-41aa-4a29-accc-bf948a7f83b5
heatmap(1:13,trange_J,otocs[12],title="H = 12",yscale=:log)

# ╔═╡ Cell order:
# ╠═df817bae-9639-44aa-bfdd-ecbc5c6b6f7e
# ╠═66e8d6f5-7bd2-4894-a9e7-823c9be8d60b
# ╠═97f67aaa-c3d4-4609-ade1-e92943b02c9c
# ╠═7af4ead4-0c87-497e-b100-566f3695ede4
# ╠═12ce7fb1-ad59-4094-ad74-c3f356e6b3be
# ╠═90bfce21-cefe-456c-8865-9f5e586f2023
# ╠═e00ef8bb-d5e3-4459-a214-c103c26a27a8
# ╠═67a4e20e-3198-4a61-be76-45c8b11265d0
# ╠═b841e087-edb1-41f0-91f8-5a9bd84367cc
# ╠═8d58599f-b79d-4160-9c8d-869554430360
# ╠═9414a48c-fe72-4c07-80fc-01d1d78c2e73
# ╠═4e203efc-bea7-47d7-ace5-50da4f43115b
# ╠═174e34f7-b948-4338-a272-5b69d65cfe50
# ╠═0a8e9c9e-185f-431b-8e8f-89eea9c17ef4
# ╠═0e9d1f52-837b-445e-9d2c-cebbb8fcdd05
# ╠═bedfd8a7-59cb-4906-9a4e-220e5a0a1a6c
# ╠═dba9941e-609c-4852-a921-1820ccb678da
# ╠═ba99f89f-527c-478e-880d-b99fd71fa67e
# ╠═0913e415-9308-4ad5-a16e-cbda2b6721a6
# ╠═98acf4a3-7a07-4039-9772-350622a3bb9a
# ╠═99cd41f2-1d32-4719-9b69-093623e5fdd4
# ╠═2fb06cd2-d0d6-46c7-bb83-544435286c51
# ╠═5c6536dd-78b9-4403-82f3-3b340916d805
# ╠═4d4ad00e-0a23-4f70-bbc1-de17c278377e
# ╠═e9554486-65cd-4a09-b2cc-1853ce075c00
# ╠═7173829a-ca92-4639-bc9c-cf545ebaa42a
# ╠═e9ca6017-41aa-4a29-accc-bf948a7f83b5
