### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ 2b39b6c0-4255-4008-8413-bda1ea656650
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 313705be-481a-11ed-17bb-f130b983a175
using LinearAlgebra,Plots,JLD2,Statistics,PlutoUI

# ╔═╡ 89e38344-a7c8-4581-ab5e-2e031d1b05a2
TableOfContents()

# ╔═╡ 01b54b4d-47dc-4a57-9490-b691de6fd997
begin
		function disorder_mean(A,n_shots)
			return mean(A[:,:,1:n_shots,:];dims=3)[:,:,1,:]
		end
		
		function disorder_std(A,n_shots)
			return std(A[:,:,1:n_shots,:];dims=3)[:,:,1,:]
		end

		function disorder_mean(A,n_shots)
			return mean(A[:,:,1:n_shots];dims=3)[:,:,1]
		end
		
		function disorder_std(A,n_shots)
			return std(A[:,:,1:n_shots];dims=3)[:,:,1]
		end

		function state_mean(A,n_states)
			return mean(A[:,:,:,1:n_states];dims=4)[:,:,:,1]
		end

		function state_std(A,n_states)
			return std(A[:,:,:,1:n_states];dims=4)[:,:,:,1]
		end
	
		#function pos_mean(A)
		#	return mean(A;dims=2)[:,1]
		#end
		
		#function reduce_by_last(A)
			#return A[:,:,1]
		#end
end

# ╔═╡ 264929c4-375e-4146-a61c-9ce0782433d9
begin
	N = 11
	path = pwd() *"/test/"
	shots = 50
	states = 50
	trange = 0:0.1:5
	T = length(trange)
end

# ╔═╡ 9c6e9d29-8041-494a-b30b-ff042cad3fe5
md"# No disorder h=0"

# ╔═╡ da94efac-baa7-4b7f-a049-554b026da248
begin
	f_df = "7153298_N11_BS.jld2"
	f_sf = "7153301_N11_BS.jld2"
	f_m = "7153309_N11_BS.jld2"
	
	jobids_df = load(path*f_df,"jobid")
	params_df = load(path*f_df,"params")
	data_df = 2*ones(T,N,shots,states)-2*load(path*f_df,"data")

	jobids_sf = load(path*f_sf,"jobid")
	params_sf = load(path*f_sf,"params")
	data_sf = 2*ones(T,N,shots,states)-2*load(path*f_sf,"data")

	jobids_m = load(path*f_m,"jobid")
	params_m = load(path*f_m,"params")
	data_m = 2*ones(T,N,shots*states)-2*load(path*f_m,"data")

	size(data_m)
end

# ╔═╡ b6b888c7-b46f-49d0-8115-e755ec272027
plot(trange,disorder_mean(data_m,2500))

# ╔═╡ b83751a1-3ab9-41f3-b592-3c794762e8dd
maximum(data_sf-data_df)

# ╔═╡ cc9e8551-6291-41fb-b0dd-36f4c2f2849d
md"# Weak disorder h=4"

# ╔═╡ 4357fec9-fec3-428e-99d4-210bf6150695


# ╔═╡ d407860a-5fc0-4928-bbb9-a86f79e29e86


# ╔═╡ d9c874bf-37b5-4a5d-b1f0-880016ee7302
md"# Strong disorder h=12"

# ╔═╡ 9e78403a-600f-48f0-bfd7-82c3893fd63b


# ╔═╡ ff2f277f-2d0e-448a-979d-4050b35a0983


# ╔═╡ Cell order:
# ╠═313705be-481a-11ed-17bb-f130b983a175
# ╠═2b39b6c0-4255-4008-8413-bda1ea656650
# ╠═89e38344-a7c8-4581-ab5e-2e031d1b05a2
# ╠═01b54b4d-47dc-4a57-9490-b691de6fd997
# ╠═264929c4-375e-4146-a61c-9ce0782433d9
# ╠═9c6e9d29-8041-494a-b30b-ff042cad3fe5
# ╠═da94efac-baa7-4b7f-a049-554b026da248
# ╠═b6b888c7-b46f-49d0-8115-e755ec272027
# ╠═b83751a1-3ab9-41f3-b592-3c794762e8dd
# ╠═cc9e8551-6291-41fb-b0dd-36f4c2f2849d
# ╠═4357fec9-fec3-428e-99d4-210bf6150695
# ╠═d407860a-5fc0-4928-bbb9-a86f79e29e86
# ╠═d9c874bf-37b5-4a5d-b1f0-880016ee7302
# ╠═9e78403a-600f-48f0-bfd7-82c3893fd63b
# ╠═ff2f277f-2d0e-448a-979d-4050b35a0983
