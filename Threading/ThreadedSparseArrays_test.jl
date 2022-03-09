### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 6bf1cce9-c1cf-44d8-ba54-d72304349b3a
import Pkg

# ╔═╡ fcdd4cc7-9921-4bf5-b109-565483928a8d
begin
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 237c5dd8-9d58-11ec-3017-6536d35b0dcf
using SparseArrays, LinearAlgebra, Plots, KrylovKit, SpinSymmetry,ThreadedSparseArrays

# ╔═╡ f9a098b9-a730-4483-b5f7-c11fb44418ea
begin
	N = 12
	H = xxz(N,6)
	ψ0 = random_state(N)
	δt = 0.1
	trange = 0:δt:5
end

# ╔═╡ ba4c5595-b598-4d03-9dbc-3eb48ade7f05
Threads.nthreads()

# ╔═╡ 7b70bd4b-843e-4e12-8d05-1b59fd4455e5
Ht = ThreadedSparseMatrixCSC(H)

# ╔═╡ 97f34c7c-c204-4232-9e9a-f42fc349e4d0
typeof(Ht)

# ╔═╡ a5268a69-8d89-4ef5-8755-8de22a7e5fce
length(nonzeros(H))/(size(H)[1]^2)*100

# ╔═╡ 96e3da27-3e6b-45f9-90e3-ea611e699412
#Not a selfadjoint operator
begin
	T = xxz(N,6)
	T[1,2] = 2.5
	Tt = ThreadedSparseMatrixCSC(T)
	T'ψ0  == T*ψ0
end

# ╔═╡ 887209d4-ef08-4ad4-b06b-b82db4d5c903
Tt'ψ0 == Tt*ψ0

# ╔═╡ 697bf26d-2d90-4a9b-8599-dd96032d9a41
sum(abs.(Tt'ψ0-Tt*ψ0))

# ╔═╡ f313499a-cb00-4c07-8c09-10295d47f423
begin #Matrix instead of vector
	Ψ = Array{ComplexF64,2}(undef,2^N,1)
	Ψ[:,1] = ψ0
end

# ╔═╡ e9815cb8-fd22-4892-894e-95fbcc14837a
Tt*Ψ == Tt'Ψ

# ╔═╡ 05c7e37c-0b85-4881-8e3e-ff4a673951de
Ht*Ψ == Ht'Ψ

# ╔═╡ 6fdaf4c3-259c-423d-9cd9-a808532459ef
length(nonzeros(H))/(2^N*2^N)*100

# ╔═╡ a66d47d7-ffe0-4b08-adc3-f68aaabbf4a6
md"## Simple multiplication"

# ╔═╡ 10602dc6-1729-412e-af86-2eff61b904e7
Threads.nthreads()

# ╔═╡ c4ec8618-1f64-4f82-b8a6-780731d03c6e
@elapsed for i in 1:1000
	H*ψ0
end

# ╔═╡ eca7a13e-710e-402f-81e6-ac3810c23c28
@elapsed for i in 1:1000
	H'ψ0
end

# ╔═╡ d62c9283-caf0-466c-8890-8fa6c7f3617a
H'ψ0 == H*ψ0

# ╔═╡ b0c0df5c-fd88-43b7-aee0-c64f1b7c82ce
@elapsed for i in 1:1000
	Ht*ψ0
end

# ╔═╡ d0fa0671-3e08-4d8e-85d0-7f47af764d6d
@elapsed for i in 1:1000
	Ht'ψ0
end

# ╔═╡ ac68258b-49f7-40bc-93d1-57414de49139
@elapsed for i in 1:1000
	Ht'Ψ
end

# ╔═╡ 379011e4-66ba-43e6-a5ee-3ecda9d01849
Ht'*ψ0 == Ht'ψ0

# ╔═╡ 3e364247-85fb-49c2-8942-bb62e4cbbfd1
Ht'ψ0 == Ht*ψ0

# ╔═╡ 3a431456-5d0e-44bf-9e1a-221f4d0f85bf
md"## Krylov Steps"

# ╔═╡ 6b762019-60a3-4502-bda5-0851c84d33af
@elapsed for i in 1:10
	krylov_step(H,0.1,ψ0)
end

# ╔═╡ ec34dd0a-3657-4277-8373-86c4b61848ed
@elapsed for i in 1:10
	krylov_step(Ht,0.1,ψ0)
end

# ╔═╡ 90532960-4dc6-4d40-9896-e6c1617c5a7b
begin
	A = single_spin_op(σz,6,N)
	B = single_spin_op(σz,12,N)
end

# ╔═╡ b0e0b49a-d0ef-4119-afc4-f63e8fa8e867
@elapsed otoc(H,A,B,2.0,ψ0)

# ╔═╡ 5364e306-3a18-47e9-88e7-6054f1f5148c
@elapsed otoc(Ht',A,B,2.0,ψ0)

# ╔═╡ bcc12ccc-b2d6-4ef7-972e-214a71719266
typeof(Ht')

# ╔═╡ 21800790-ad55-405d-8890-fcac9b1250bd
md"# Large Disorder"

# ╔═╡ 8b79d997-493f-4113-b749-9a41cf4f6e19
begin
	L = 17
	Al = single_spin_op(σz,6,L)
	Bl = single_spin_op(σz,12,L)
	Hl = xxz(L,6.0) + field_term(12.0,L)
	Hlt = ThreadedSparseMatrixCSC(Hl)'
	ψl = random_state(L)
end

# ╔═╡ 3bef4b91-9356-4928-a479-5913e748974a
@elapsed otoc(Hlt,Al,Bl,1.0,ψl)

# ╔═╡ 13ce5125-aff3-4dca-828b-8df860c40c50
begin
	hs = 0.0:2.0:14.0
	time = Vector{Float64}(undef,length(hs))
	time2 = Vector{Float64}(undef,length(hs))
	time3 = Vector{Float64}(undef,length(hs))
	time4 = Vector{Float64}(undef,length(hs))
	for (i,h) in enumerate(hs)
		Hl2 = xxz(L,6.0) + field_term(h,L)
		Hlt2 = ThreadedSparseMatrixCSC(Hl2)'
		time[i] = @elapsed otoc(Hl2,Al,Bl,1.0,ψl)
		print("1 done\n")
		time2[i] = @elapsed otoc(Hlt2,Al,Bl,1.0,ψl)
		print("1 done\n")
		time3[i] = @elapsed otoc(Hl2,Al,Bl,5.0,ψl)
		print("1 done\n")
		time4[i] = @elapsed otoc(Hlt2,Al,Bl,5.0,ψl)
		print("1 done\n")
	end
end

# ╔═╡ 96303cc0-5aba-41e8-93f2-76a056d934e7
time

# ╔═╡ 12d10c37-9f73-49bc-8d4b-c06bc4c8daaf
54.8897/16.5367 #17

# ╔═╡ 8c1f9627-a17f-4e7d-8339-c61b415574e0
13.5212/6.59397 #16

# ╔═╡ b6e509c9-95fa-43b5-aae8-0f576b39e891
1.93161/0.768888 #14

# ╔═╡ 04c3c1a2-3d24-49df-b349-a2bb4f7b792f
0.176421/0.0888555 #12

# ╔═╡ 987b400c-97bb-4a18-a346-28843f5d9721
begin
	plot(hs,time,label="normal, t=1",xlabel="h",ylabel="t")
	plot!(hs,time2,label="threaded, t=1")
end

# ╔═╡ 01f644a5-85c5-4dfe-bbc3-331790780c75
time3

# ╔═╡ de44342e-dd19-4e7e-a5a9-85e4fa033e6a
begin
		plot(hs,time3,label="normal, t=5",xlabel="h",ylabel="t")
		plot!(hs,time4,label="threaded, t=5")
end

# ╔═╡ Cell order:
# ╠═6bf1cce9-c1cf-44d8-ba54-d72304349b3a
# ╠═fcdd4cc7-9921-4bf5-b109-565483928a8d
# ╠═237c5dd8-9d58-11ec-3017-6536d35b0dcf
# ╠═f9a098b9-a730-4483-b5f7-c11fb44418ea
# ╠═ba4c5595-b598-4d03-9dbc-3eb48ade7f05
# ╠═7b70bd4b-843e-4e12-8d05-1b59fd4455e5
# ╠═97f34c7c-c204-4232-9e9a-f42fc349e4d0
# ╠═a5268a69-8d89-4ef5-8755-8de22a7e5fce
# ╠═96e3da27-3e6b-45f9-90e3-ea611e699412
# ╠═887209d4-ef08-4ad4-b06b-b82db4d5c903
# ╠═697bf26d-2d90-4a9b-8599-dd96032d9a41
# ╠═f313499a-cb00-4c07-8c09-10295d47f423
# ╠═e9815cb8-fd22-4892-894e-95fbcc14837a
# ╠═05c7e37c-0b85-4881-8e3e-ff4a673951de
# ╠═6fdaf4c3-259c-423d-9cd9-a808532459ef
# ╠═a66d47d7-ffe0-4b08-adc3-f68aaabbf4a6
# ╠═10602dc6-1729-412e-af86-2eff61b904e7
# ╠═c4ec8618-1f64-4f82-b8a6-780731d03c6e
# ╠═eca7a13e-710e-402f-81e6-ac3810c23c28
# ╠═d62c9283-caf0-466c-8890-8fa6c7f3617a
# ╠═b0c0df5c-fd88-43b7-aee0-c64f1b7c82ce
# ╠═d0fa0671-3e08-4d8e-85d0-7f47af764d6d
# ╠═ac68258b-49f7-40bc-93d1-57414de49139
# ╠═379011e4-66ba-43e6-a5ee-3ecda9d01849
# ╠═3e364247-85fb-49c2-8942-bb62e4cbbfd1
# ╠═3a431456-5d0e-44bf-9e1a-221f4d0f85bf
# ╠═6b762019-60a3-4502-bda5-0851c84d33af
# ╠═ec34dd0a-3657-4277-8373-86c4b61848ed
# ╠═90532960-4dc6-4d40-9896-e6c1617c5a7b
# ╠═b0e0b49a-d0ef-4119-afc4-f63e8fa8e867
# ╠═5364e306-3a18-47e9-88e7-6054f1f5148c
# ╠═bcc12ccc-b2d6-4ef7-972e-214a71719266
# ╠═21800790-ad55-405d-8890-fcac9b1250bd
# ╠═8b79d997-493f-4113-b749-9a41cf4f6e19
# ╠═3bef4b91-9356-4928-a479-5913e748974a
# ╠═13ce5125-aff3-4dca-828b-8df860c40c50
# ╠═96303cc0-5aba-41e8-93f2-76a056d934e7
# ╠═12d10c37-9f73-49bc-8d4b-c06bc4c8daaf
# ╠═8c1f9627-a17f-4e7d-8339-c61b415574e0
# ╠═b6e509c9-95fa-43b5-aae8-0f576b39e891
# ╠═04c3c1a2-3d24-49df-b349-a2bb4f7b792f
# ╠═987b400c-97bb-4a18-a346-28843f5d9721
# ╠═01f644a5-85c5-4dfe-bbc3-331790780c75
# ╠═de44342e-dd19-4e7e-a5a9-85e4fa033e6a
