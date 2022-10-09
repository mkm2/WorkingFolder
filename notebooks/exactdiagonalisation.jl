### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ c8339c8e-383b-11ed-20f5-930afe4d75a8
import Pkg

# ╔═╡ 36cad27b-2ebb-4ca1-abc0-8e4b21595e4a
begin
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 843165f8-c92b-4ff8-ace1-b20c1f9bf891
using SparseArrays, LinearAlgebra, Plots, SpinSymmetry, BenchmarkTools

# ╔═╡ b2c2cedc-c382-450b-84d8-026f3c1a543e
using Tullio, PlutoUI #,MatrixChainMultiply

# ╔═╡ d6b82f30-3e38-492d-88e9-8bdd414935fd
TableOfContents()

# ╔═╡ ed61ced6-61fe-4d37-bca1-ee633abee735
html"""<style>main {max-width: 60%;}</style>"""

# ╔═╡ 484e4df3-aebb-4ea7-bc1a-61262d252872
begin
	N = 9
	i = div(N,2)+1
	k = div(N-1,2)+1 #largest sector
	d = basissize(symmetrized_basis(N,k))
	H = xxz(N,6)
	δt = 0.1
	trange = 0.1:δt:5
	A = single_spin_op(σz,i,N)
end

# ╔═╡ afd40ad1-4679-4814-aa40-d56537065f33
md"# Exact diagonalisation"

# ╔═╡ 75be913e-2fdd-4b8d-9611-6aca259ac37c
evals, evecs = eigen(Hermitian(Matrix(H)))

# ╔═╡ 4c9aac44-5422-4077-a204-c715cc230402
begin
	Aeig = evecs'*A*evecs
	B = single_spin_op(σz,4,N)
	Beig = evecs'*B*evecs
	em = exp(-im*Diagonal(evals))
	ψ0 = random_state(N)
	ψeig = evecs'*ψ0
end

# ╔═╡ 26dd5aae-561f-4020-9aa8-a5c8bd652aad
function unit(size,index)
	vec = zeros(Float64,size)
	vec[index] = 1.0
	return vec
end

# ╔═╡ 4e988e7f-e445-4ab3-843b-ebb4273d8b47
begin
	function Fψ(t::Float64,A::Matrix{Float64},B::Matrix{Float64},evals::Vector{Float64},ψ)
		d = size(A)[1]
		eigmt = exp(-im*Diagonal(evals)*t)
		return ψ'*eigmt'*A*eigmt*B*eigmt'*A*eigmt*B*ψ #much faster than trace!
	end

	function Frand(t,A,B,evals,N)
		return Fψ(t,A,B,evals,random_state(N))
	end

	Frandmean(t,A,B,evals,N,s) = mean(Frand(t,A,B,evals,N) for i in 1:s)

	function Ftr(t::Float64,A::Matrix{Float64},B::Matrix{Float64},evals::Vector{Float64})
		d = size(A)[1]
		eigmt = exp(-im*Diagonal(evals)*t)
		return tr(eigmt'*A*eigmt*B*eigmt'*A*eigmt*B)
	end
	
	function Fsum(t::Float64,A::Matrix{Float64},B::Matrix{Float64},evals::Vector{Float64})
		d = size(A)[1]
		eigmt = exp(-im*Diagonal(evals)*t)
		tmp = eigmt'*A*eigmt*B*eigmt'*A*eigmt*B
		f = 0
		for i in 1:d
			f += unit(d,i)'tmp*unit(d,i)
		end
		return f
	end

	function Fmatchain(t::Float64,A::Matrix{Float64},B::Matrix{Float64},evals::Vector{Float64})
		d = size(A)[1]
		eigmt = exp(-im*Diagonal(evals)*t)
		return tr(matrixchainmultiply(eigmt',A,eigmt,B,eigmt',A,eigmt,B))
	end

	function FTullio(t::Float64,A::Matrix{Float64},B::Matrix{Float64},evals::Vector{Float64})
		d = size(A)[1]
		eigmt = exp(-im*Diagonal(evals)*t) 
		return @tullio f:= A[n,k] * B[k,m] * A[m,l] * B[l,n] * exp(im*(evals[n]-evals[k])*t) * exp(im*(evals[m]-evals[l])*t)
	end
end

# ╔═╡ 3e7f2f54-fb02-40fa-b635-e8ae951b5641
@btime otoc(H,A,B,1.2,ψ0)

# ╔═╡ 8f31b022-2b10-439e-9ef7-80df1588965a
@benchmark Ftr(1.0,Aeig,Beig,evals)/2^N

# ╔═╡ 41a1d686-38a6-4a71-9437-5c6a3a549be5
@benchmark Fsum(1.0,Aeig,Beig,evals)/2^N

# ╔═╡ 823d0291-d192-44f0-9e73-8d1e2e9d7df7
@benchmark Fmatchain(1.0,Aeig,Beig,evals)/2^N

# ╔═╡ c42063b6-9b00-4aed-9887-d71359cd1f52
#FTullio(1.0,Aeig,Beig,evals)/2^N

# ╔═╡ b485f8fa-659d-439e-ab4a-e774d082c497
@benchmark Frand(1.0,Aeig,Beig,evals,N)

# ╔═╡ f5597b69-964a-45f7-b112-a88eab18fc53
@benchmark Frandmean(1.0,Aeig,Beig,evals,N,10)

# ╔═╡ c7ac330c-d6e3-489c-9426-4973a14921d0
begin #B not in eigenbasis
	function FψQ(t::Float64,A::Matrix{Float64},B::AbstractArray{Float64},Q,evals::Vector{Float64},ψ)
			d = size(A)[1]
			eigmt = exp(-im*Diagonal(evals)*t)
			return ψ'*eigmt'*A*eigmt*Q'*B*Q*eigmt'*A*eigmt*Q'*B*Q*ψ #much faster than trace!
		end
	FψQ(t::Float64,A::Matrix{Float64},B::AbstractArray{Float64},evals::Vector{Float64},ψ,l) = FψQ(t::Float64,A::Matrix{Float64},B::AbstractArray{Float64},Diagonal(ones(2^l)),evals::Vector{Float64},ψ)
end

# ╔═╡ d677aed2-debf-4578-bac3-97902682e06f
@benchmark FψQ(1.0,Aeig,B,evecs,evals,ψ0)

# ╔═╡ 6029cc6c-a61e-4b6d-b60b-dd9f639fc8ca
@benchmark FψQ(1.0,Aeig,Beig,evals,ψ0,N)

# ╔═╡ 4cdd0c2e-3ea6-4a2a-b282-76bc74922f89
begin
	S = 100
	test = zeros(S)
	ref = Ftr(1.0,Aeig,Beig,evals)/2^N
	for i in 1:S
		test[i] = abs(Frandmean(1.0,Aeig,Beig,evals,N,i)-ref)/abs(ref)
	end
end

# ╔═╡ 40b45ebd-1ccc-4449-ac4e-c2397fa14c75
plot(1:S,test*100)

# ╔═╡ 417f8665-de9f-4276-8752-6e878b530e74
begin
	Atest = convert(Matrix{ComplexF64},Aeig)
	Btest = convert(Matrix{ComplexF64},B)
	
end

# ╔═╡ ac0dda2b-f27e-47a5-a312-0be44d2d722d
begin
	corr = zeros(length(trange))
	corr2 = similar(corr)
	corr3 = similar(corr)
	for (i,t) in enumerate(trange)
		corr[i] = 2*(1-otoc(H,A,B,t,ψ0))
		corr2[i] = 2*(1-real(FψQ(t,Aeig,Beig,evals,ψeig,N)))
		corr3[i] = 2*(1-otoc_edψ(Atest,Btest,evals,evecs,t,ψeig))
	end
end

# ╔═╡ 7eb5529e-fcb1-468d-8b0d-ad34d8df807c
begin
	plot(trange,corr)
	plot!(trange,corr2)
	plot!(trange,corr3)
end

# ╔═╡ 5f75543f-f484-41b9-85f7-839a30841725
otoc_edψ(Atest,Btest,evals,evecs,1.0,ψ0)

# ╔═╡ 1cb0bb29-52c9-4daf-adbf-94e0432469e2
convert(Vector{Float64},evals)

# ╔═╡ 67134d83-709b-4147-987b-32754e2c162c
Diagonal(ones(2^N)) isa AbstractArray

# ╔═╡ 2d7998f4-912e-4ce9-b2ca-9080266701aa
md"# Speed of MatrixChainMultiply"

# ╔═╡ 76666161-55e9-4560-9b5b-5d8a36ea1a6e
begin
	a = rand(1000,1000)
	b = rand(1000,100)
	c = rand(100, 500)
	D = rand(500)
end

# ╔═╡ 89676e1a-4fd4-435e-a7f8-1b1223ec4396
@benchmark result1 = matrixchainmultiply(a,b,c,D)

# ╔═╡ fa40cde6-96a5-4fb1-923d-a961ea0b417d
@benchmark result2 = *(a,b,c,D)

# ╔═╡ Cell order:
# ╠═c8339c8e-383b-11ed-20f5-930afe4d75a8
# ╠═36cad27b-2ebb-4ca1-abc0-8e4b21595e4a
# ╠═843165f8-c92b-4ff8-ace1-b20c1f9bf891
# ╠═b2c2cedc-c382-450b-84d8-026f3c1a543e
# ╠═d6b82f30-3e38-492d-88e9-8bdd414935fd
# ╠═ed61ced6-61fe-4d37-bca1-ee633abee735
# ╠═484e4df3-aebb-4ea7-bc1a-61262d252872
# ╠═afd40ad1-4679-4814-aa40-d56537065f33
# ╠═75be913e-2fdd-4b8d-9611-6aca259ac37c
# ╠═4c9aac44-5422-4077-a204-c715cc230402
# ╠═4e988e7f-e445-4ab3-843b-ebb4273d8b47
# ╠═26dd5aae-561f-4020-9aa8-a5c8bd652aad
# ╠═3e7f2f54-fb02-40fa-b635-e8ae951b5641
# ╠═8f31b022-2b10-439e-9ef7-80df1588965a
# ╠═41a1d686-38a6-4a71-9437-5c6a3a549be5
# ╠═823d0291-d192-44f0-9e73-8d1e2e9d7df7
# ╠═c42063b6-9b00-4aed-9887-d71359cd1f52
# ╠═b485f8fa-659d-439e-ab4a-e774d082c497
# ╠═f5597b69-964a-45f7-b112-a88eab18fc53
# ╠═c7ac330c-d6e3-489c-9426-4973a14921d0
# ╠═d677aed2-debf-4578-bac3-97902682e06f
# ╠═6029cc6c-a61e-4b6d-b60b-dd9f639fc8ca
# ╠═4cdd0c2e-3ea6-4a2a-b282-76bc74922f89
# ╠═40b45ebd-1ccc-4449-ac4e-c2397fa14c75
# ╠═ac0dda2b-f27e-47a5-a312-0be44d2d722d
# ╠═7eb5529e-fcb1-468d-8b0d-ad34d8df807c
# ╠═417f8665-de9f-4276-8752-6e878b530e74
# ╠═5f75543f-f484-41b9-85f7-839a30841725
# ╠═1cb0bb29-52c9-4daf-adbf-94e0432469e2
# ╠═67134d83-709b-4147-987b-32754e2c162c
# ╠═2d7998f4-912e-4ce9-b2ca-9080266701aa
# ╠═76666161-55e9-4560-9b5b-5d8a36ea1a6e
# ╠═89676e1a-4fd4-435e-a7f8-1b1223ec4396
# ╠═fa40cde6-96a5-4fb1-923d-a961ea0b417d
