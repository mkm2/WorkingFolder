### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 96212858-0439-11ed-23f0-377150c3ddaa
import Pkg

# ╔═╡ 60775287-1d71-43ca-9a30-8e9d44e81a84
begin
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ e93c4236-7a19-41a9-be62-e5feef4f4993
using SparseArrays, LinearAlgebra, Plots, KrylovKit, SpinSymmetry

# ╔═╡ ff37cea5-6bef-4084-9d96-0dfd9ec40f1b
begin
	N = 12
	H = xxz(N,6)
	ψ0 = normalize!(ones(2^N))
	tmax = 1
	trange = 0:0.1:5
end

# ╔═╡ 633277f4-3cee-41de-9ad1-ff8123d8bf12
symmetrize_operator(H,N,div(N-1,2)+1)

# ╔═╡ 8c9465d2-a548-4f2e-8114-75976479007b
for i in 1:6
	print(i,div(i-1,2)+1)
end #correct one

# ╔═╡ 1b458a9f-2b9d-4f40-861d-e2d0ca83d679
begin
	function random_state2(N::Int, d::Int)
	    return normalize!(randn(ComplexF64,d))
	end
	random_state2(N::Int) = random_state2(N, 2^N)
end

# ╔═╡ 5fdec16a-0480-4af1-b41a-7cbef0366292
begin
	function random_bitstring_state2(N::Int, d::Int)
	    x = zeros(ComplexF64,d)
	    ind = rand(1:d)
	    x[ind] = 1
	    return x
	end
	random_bitstring_state2(N::Int) = random_bitstring_state2(N, 2^N)
end

# ╔═╡ e8b823f9-a3bb-44c1-b511-667a6bce7f42
begin
	x = random_product_state(3)
	x
end

# ╔═╡ 8bc6e8cb-3d15-4834-85ff-b81d40146074
normalize!(symmetrize_state(x,3,2))

# ╔═╡ 64112ad0-e34e-4c4a-a1a5-7aeb52044826
begin
	k = div(N-1,2)+1
	d = basissize(symmetrized_basis(N,k))
	H_eff = symmetrize_operator(H,N,k)
	ψ_eff = random_state2(N,d)
	op1 = single_spin_op(σz,5,N)
	op1_eff = symmetrize_operator(op1,N,k)
	op2 = single_spin_op(σz,1,N)
	op2_eff = symmetrize_operator(op2,N,k)
end

# ╔═╡ aa2ffde7-f93d-4205-be09-6c834f61d884
begin
	corr = zeros(Float64,length(trange))
	@elapsed for (ti,t) in enumerate(trange)
	    corr[ti] = 2-2*otoc(H_eff, op1_eff, op2_eff, t, random_state2(N,d))
	end
end

# ╔═╡ 5f6f0fae-7a02-4c2f-8e89-a69de24487fd
plot(trange,corr)

# ╔═╡ Cell order:
# ╠═96212858-0439-11ed-23f0-377150c3ddaa
# ╠═60775287-1d71-43ca-9a30-8e9d44e81a84
# ╠═e93c4236-7a19-41a9-be62-e5feef4f4993
# ╠═ff37cea5-6bef-4084-9d96-0dfd9ec40f1b
# ╠═633277f4-3cee-41de-9ad1-ff8123d8bf12
# ╠═8c9465d2-a548-4f2e-8114-75976479007b
# ╠═1b458a9f-2b9d-4f40-861d-e2d0ca83d679
# ╠═5fdec16a-0480-4af1-b41a-7cbef0366292
# ╠═e8b823f9-a3bb-44c1-b511-667a6bce7f42
# ╠═8bc6e8cb-3d15-4834-85ff-b81d40146074
# ╠═64112ad0-e34e-4c4a-a1a5-7aeb52044826
# ╠═aa2ffde7-f93d-4205-be09-6c834f61d884
# ╠═5f6f0fae-7a02-4c2f-8e89-a69de24487fd
