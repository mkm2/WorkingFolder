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
		function disorder_mean(A::Array{Float64,4},n_shots)
			return mean(A[:,:,1:n_shots,:];dims=3)[:,:,1,:]
		end
		
		function disorder_std(A::Array{Float64,4},n_shots)
			return std(A[:,:,1:n_shots,:];dims=3)[:,:,1,:]
		end

		function disorder_mean(A::Array{Float64,3},n_shots)
			return mean(A[:,:,1:n_shots];dims=3)[:,:,1]
		end
		
		function disorder_std(A::Array{Float64,3},n_shots)
			return std(A[:,:,1:n_shots];dims=3)[:,:,1]
		end

		function state_mean(A::Array{Float64,4},n_states)
			return mean(A[:,:,:,1:n_states];dims=4)[:,:,:,1]
		end

		function state_std(A::Array{Float64,4},n_states)
			return std(A[:,:,:,1:n_states];dims=4)[:,:,:,1]
		end

		function state_mean(A::Array{Float64,3},n_states)
			return mean(A[:,:,1:n_states];dims=3)[:,:,1]
		end

		function state_std(A::Array{Float64,3},n_states)
			return std(A[:,:,1:n_states];dims=3)[:,:,1]
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
	N = 13
	path = pwd() *"/test/"
	shots = 50
	states = 50
	trange = 0:0.1:5
	T = length(trange)
end

# ╔═╡ 421be086-0122-4a72-98b3-c3914f032c71
begin
	ψs = Vector{Vector{ComplexF64}}(undef,10)
		for ind in 1:10
			ψs[ind] = random_state(N)
		end
end

# ╔═╡ dc6e8621-b575-4014-987b-ae9f38431c29
ψs

# ╔═╡ dbd00835-2ca3-4dc9-875b-71fb2dcc03d1
N

# ╔═╡ 9c6e9d29-8041-494a-b30b-ff042cad3fe5
md"# No disorder h=0"

# ╔═╡ da94efac-baa7-4b7f-a049-554b026da248
begin
	f_df = "7301823_N13_BS.jld2"
	f_sf = "7273738_N13_BS.jld2"
	f_m = "7273742_N13_BS.jld2"
	f_ref = "7296575_N13_RS.jld2"
	
	jobids_df = load(path*f_df,"jobid")
	params_df = load(path*f_df,"params")
	data_df = 2*ones(T,N,shots,states)-2*load(path*f_df,"data")

	jobids_sf = load(path*f_sf,"jobid")
	params_sf = load(path*f_sf,"params")
	data_sf = 2*ones(T,N,shots,states)-2*load(path*f_sf,"data")

	jobids_m = load(path*f_m,"jobid")
	params_m = load(path*f_m,"params")
	data_m = 2*ones(T,N,shots*states)-2*load(path*f_m,"data")

	jobids_ref = load(path*f_ref,"jobid")
	params_ref = load(path*f_ref,"params")
	data_ref = 2*ones(T,N,100,10)-2*load(path*f_ref,"data")

	size(data_ref)
end

# ╔═╡ 56b08a04-4ae1-43b6-ae5a-32e997357f3a
ref = disorder_mean(state_mean(data_ref,10),100)

# ╔═╡ f8b8c837-ed9e-43f6-a2fb-b1e6908c010b
begin
	sample_states = 5
	sample_disorder = 10
end

# ╔═╡ d37571e5-9e66-4410-815f-71227cf15598
begin
	df_mean_byd = Vector{Matrix{Float64}}(undef,shots)
	sf_mean_bys = Vector{Matrix{Float64}}(undef,states)
	m_mean_bysample = Vector{Matrix{Float64}}(undef,shots*states)
	for i in 1:shots
		df_mean_byd[i] = disorder_mean(state_mean(data_df,sample_states),i)
	end
	for i in 1:states
		sf_mean_bys[i] = state_mean(disorder_mean(data_sf,sample_disorder),i)
	end
	for i in 1:shots*states
		m_mean_bysample[i] = disorder_mean(data_m,i)
	end
end

# ╔═╡ 879c6440-2888-4770-93ff-7e33fe8cf1ee
begin
	df_meandeviation = Vector{Float64}(undef,shots)
	sf_meandeviation = Vector{Float64}(undef,states)
	m_meandeviation = Vector{Float64}(undef,shots*states)
	
	for i in 1:shots
		df_meandeviation[i] = sum(abs.(df_mean_byd[i]-ref))
	end
	for i in 1:states
		sf_meandeviation[i] = sum(abs.(sf_mean_bys[i]-ref))
	end
	for i in 1:shots*states
		m_meandeviation[i] = sum(abs.(m_mean_bysample[i]-ref))
	end
end

# ╔═╡ f1b5be26-ef2d-4b1a-8284-797171679664
m_meandeviation/(51*13)

# ╔═╡ d83ee2b2-4163-43f1-b0ee-3bcc39a7b580
sample_states * (1:50)

# ╔═╡ bb450ebe-b836-4143-be7a-dccf930aad4c


# ╔═╡ fc2c15cd-156e-4c89-a90d-330d88d1cd25
begin
	plot(sample_states * (1:50),df_meandeviation,yaxis=:log)
	plot!(sample_disorder * (1:50),sf_meandeviation)
	plot!(1:2500,m_meandeviation)
end

# ╔═╡ cc9e8551-6291-41fb-b0dd-36f4c2f2849d
md"# Weak disorder h=4"

# ╔═╡ d9c874bf-37b5-4a5d-b1f0-880016ee7302
md"# Strong disorder h=12"

# ╔═╡ ff2f277f-2d0e-448a-979d-4050b35a0983


# ╔═╡ 189594bf-20cc-4dd5-9583-cbe930a04c4a
begin
	df_meandeviation_weak = Vector{Float64}(undef,shots)
	sf_meandeviation_weak = Vector{Float64}(undef,states)
	m_meandeviation_weak = Vector{Float64}(undef,shots*states)
	
	for i in 1:shots
		df_meandeviation_weak[i] = sum(abs.(df_mean_byd_weak[i]-ref))
	end
	for i in 1:states
		sf_meandeviation_weak[i] = sum(abs.(sf_mean_bys_weak[i]-ref))
	end
	for i in 1:shots*states
		m_meandeviation_weak[i] = sum(abs.(m_mean_bysample_weak[i]-ref))
	end
end

# ╔═╡ cf45a056-dac7-4f2b-b972-05764077f107
begin
	f_df_weak = "XXX_N13_BS.jld2"
	f_sf_weak = "XXX_N13_BS.jld2"
	f_m_weak = "XXX_N13_BS.jld2"
	f_ref_weak = "XXX_N13_RS.jld2"
	
	jobids_df_weak = load(path*f_df_weak,"jobid")
	params_df_weak = load(path*f_df_weak,"params")
	data_df_weak = 2*ones(T,N,shots,states)-2*load(path*f_df_weak,"data")

	jobids_sf_weak = load(path*f_sf_weak,"jobid")
	params_sf_weak = load(path*f_sf_weak,"params")
	data_sf_weak = 2*ones(T,N,shots,states)-2*load(path*f_sf_weak,"data")

	jobids_m_weak = load(path*f_m_weak,"jobid")
	params_m_weak = load(path*f_m_weak,"params")
	data_m_weak = 2*ones(T,N,shots*states)-2*load(path*f_m_weak,"data")

	jobids_ref_weak = load(path*f_ref_weak,"jobid")
	params_ref_weak = load(path*f_ref_weak,"params")
	data_ref_weak = 2*ones(T,N,100,10)-2*load(path*f_ref_weak,"data")

	size(data_ref_weak)
end

# ╔═╡ 4357fec9-fec3-428e-99d4-210bf6150695
begin
	df_mean_byd_weak = Vector{Matrix{Float64}}(undef,shots)
	sf_mean_bys_weak = Vector{Matrix{Float64}}(undef,states)
	m_mean_bysample_weak = Vector{Matrix{Float64}}(undef,shots*states)
	for i in 1:shots
		df_mean_byd_weak[i] = disorder_mean(state_mean(data_df_weak,sample_states),i)
	end
	for i in 1:states
		sf_mean_bys_weak[i] = state_mean(disorder_mean(data_sf_weak,sample_disorder),i)
	end
	for i in 1:shots*states
		m_mean_bysample_weak[i] = disorder_mean(data_m_weak,i)
	end
end

# ╔═╡ 0c598c0f-ac2c-41c8-8174-c968a2cb6540
begin
	f_df_weak = "XXX_N13_BS.jld2"
	f_sf_weak = "XXX_N13_BS.jld2"
	f_m_weak = "XXX_N13_BS.jld2"
	f_ref_weak = "XXX_N13_RS.jld2"
	
	jobids_df_weak = load(path*f_df_weak,"jobid")
	params_df_weak = load(path*f_df_weak,"params")
	data_df_weak = 2*ones(T,N,shots,states)-2*load(path*f_df_weak,"data")

	jobids_sf_weak = load(path*f_sf_weak,"jobid")
	params_sf_weak = load(path*f_sf_weak,"params")
	data_sf_weak = 2*ones(T,N,shots,states)-2*load(path*f_sf_weak,"data")

	jobids_m_weak = load(path*f_m_weak,"jobid")
	params_m_weak = load(path*f_m_weak,"params")
	data_m_weak = 2*ones(T,N,shots*states)-2*load(path*f_m_weak,"data")

	jobids_ref_weak = load(path*f_ref_weak,"jobid")
	params_ref_weak = load(path*f_ref_weak,"params")
	data_ref_weak = 2*ones(T,N,100,10)-2*load(path*f_ref_weak,"data")

	size(data_ref_weak)
end

# ╔═╡ d407860a-5fc0-4928-bbb9-a86f79e29e86
begin
	df_meandeviation_weak = Vector{Float64}(undef,shots)
	sf_meandeviation_weak = Vector{Float64}(undef,states)
	m_meandeviation_weak = Vector{Float64}(undef,shots*states)
	
	for i in 1:shots
		df_meandeviation_weak[i] = sum(abs.(df_mean_byd_weak[i]-ref))
	end
	for i in 1:states
		sf_meandeviation_weak[i] = sum(abs.(sf_mean_bys_weak[i]-ref))
	end
	for i in 1:shots*states
		m_meandeviation_weak[i] = sum(abs.(m_mean_bysample_weak[i]-ref))
	end
end

# ╔═╡ b47b5869-fba8-458b-9b8a-03ba89dcef9b
begin
	df_mean_byd_weak = Vector{Matrix{Float64}}(undef,shots)
	sf_mean_bys_weak = Vector{Matrix{Float64}}(undef,states)
	m_mean_bysample_weak = Vector{Matrix{Float64}}(undef,shots*states)
	for i in 1:shots
		df_mean_byd_weak[i] = disorder_mean(state_mean(data_df_weak,sample_states),i)
	end
	for i in 1:states
		sf_mean_bys_weak[i] = state_mean(disorder_mean(data_sf_weak,sample_disorder),i)
	end
	for i in 1:shots*states
		m_mean_bysample_weak[i] = disorder_mean(data_m_weak,i)
	end
end

# ╔═╡ Cell order:
# ╠═313705be-481a-11ed-17bb-f130b983a175
# ╠═2b39b6c0-4255-4008-8413-bda1ea656650
# ╠═89e38344-a7c8-4581-ab5e-2e031d1b05a2
# ╠═421be086-0122-4a72-98b3-c3914f032c71
# ╠═dbd00835-2ca3-4dc9-875b-71fb2dcc03d1
# ╠═dc6e8621-b575-4014-987b-ae9f38431c29
# ╠═01b54b4d-47dc-4a57-9490-b691de6fd997
# ╠═264929c4-375e-4146-a61c-9ce0782433d9
# ╠═9c6e9d29-8041-494a-b30b-ff042cad3fe5
# ╠═da94efac-baa7-4b7f-a049-554b026da248
# ╠═56b08a04-4ae1-43b6-ae5a-32e997357f3a
# ╠═f8b8c837-ed9e-43f6-a2fb-b1e6908c010b
# ╠═d37571e5-9e66-4410-815f-71227cf15598
# ╠═879c6440-2888-4770-93ff-7e33fe8cf1ee
# ╠═f1b5be26-ef2d-4b1a-8284-797171679664
# ╠═d83ee2b2-4163-43f1-b0ee-3bcc39a7b580
# ╠═bb450ebe-b836-4143-be7a-dccf930aad4c
# ╠═fc2c15cd-156e-4c89-a90d-330d88d1cd25
# ╠═cc9e8551-6291-41fb-b0dd-36f4c2f2849d
# ╠═cf45a056-dac7-4f2b-b972-05764077f107
# ╠═4357fec9-fec3-428e-99d4-210bf6150695
# ╠═d407860a-5fc0-4928-bbb9-a86f79e29e86
# ╠═d9c874bf-37b5-4a5d-b1f0-880016ee7302
# ╠═0c598c0f-ac2c-41c8-8174-c968a2cb6540
# ╠═b47b5869-fba8-458b-9b8a-03ba89dcef9b
# ╠═189594bf-20cc-4dd5-9583-cbe930a04c4a
# ╠═ff2f277f-2d0e-448a-979d-4050b35a0983
