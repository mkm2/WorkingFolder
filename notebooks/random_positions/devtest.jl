### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ e0bdd126-9e5f-11ec-1174-e7a679707f91
import Pkg

# ╔═╡ 452127a7-c67f-4edd-a550-b87ee1392f5b
begin
	Pkg.add(url="https://github.com/abraemer/XXZNumerics.jl")
	using XXZNumerics
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ e4052cfb-5950-4b23-a6c8-0e81ee531fa7
using SparseArrays, LinearAlgebra, Plots, KrylovKit, SpinSymmetry,ThreadedSparseArrays

# ╔═╡ bf8244ec-4cab-47a0-9481-15974ed7fd96
using PlutoUI

# ╔═╡ 3baa2897-8c0a-48bf-a36d-e99b70723c43
const Maybe{T} = Union{Missing, T} where T

# ╔═╡ 74e362cd-fca3-4e54-9e3c-06f3e6f8ea3a
md"# Geometry and Positions"

# ╔═╡ effee07f-d0fa-4bee-a741-62ac8339fefe
begin #from SimLib, no changes made
	const GEOMETRIES = [:box, :box_pbc, :noisy_chain, :noisy_chain_pbc]
	
	# The OBC Volume coeffs should probably be a bit bigger to accommodate for the extra space outside
	# Right now results won't really be comparable, but we will not need this (right now that is)
	const VOLUME_COEFFS = (; sphere=[2.0, π, 4/3*π], box=[1.0,1.0,1.0], box_pbc=[1.0,1.0,1.0])
	
	parse_geometry(s::AbstractString) = Symbol(lowercase(s))
	parse_geometry(s::Symbol) = s
	
	## General idea:
	## set r_bl=1 and scale volume via the density rho
	## rho is defined as the ratio (blockaded volume) / (total volume)
	## total volume = VOLUME_FACTOR[geom][d]*L**d
	## blockaded volume = N*VOLUME_FACTOR[sphere][d]*1**d
	length_from_density(geometry, ρ, N, dim) = (1/ρ * N * VOLUME_COEFFS.sphere[dim] / VOLUME_COEFFS[geometry][dim])^(1/dim)
	
	geometry_from_density(geom::Symbol, ρ, N, dim) = geometry_from_density(Val(geom), ρ, N, dim)
	geometry_from_density(::Val{:box}, ρ, N, dim) = Box(ones(dim)*length_from_density(:box, ρ, N, dim))
	geometry_from_density(::Val{:box_pbc}, ρ, N, dim) = BoxPBC(ones(dim)*length_from_density(:box_pbc, ρ, N, dim))
	function geometry_from_density(::Val{:noisy_chain}, ρ, N, dim)
	    dim == 1 || error("Chain is only 1D!")
	    # know that rho = 2/spacing, since r_bl = 1, and V = N*spacing
	    spacing = 2/ρ
	    σ = 1.5*(spacing-1) # this ensures ~90% success rate when generating positions -> see demo/rho_scaling.ipynb
	    NoisyChain(N, spacing, σ)
	end
	
	function geometry_from_density(::Val{:noisy_chain_pbc}, ρ, N, dim)
	    dim == 1 || error("Chain is only 1D!")
	    # know that rho = 2/spacing, since r_bl = 1, and V = N*spacing
	    spacing = 2/ρ
	    σ = 1.5*(spacing-1) # this ensures ~90% success rate when generating positions -> see demo/rho_scaling.ipynb
	    NoisyChainPBC(N, spacing, σ)
	end
	
	# Determined by testing - within these ranges the position sampling converges almost certainly
	const SAFE_RHO_RANGES = (;box     = [(0.1, 1.25), (0.1, 1.7), (0.1,2.05)],
	                          box_pbc = [(0.1, 1.25), (0.1, 1.7), (0.1,2.05)],
	                          noisy_chain_pbc = [(1.0, 1.99)],
	                          noisy_chain     = [(1.0, 1.99)])
end

# ╔═╡ 66cccf86-5d37-4e38-93d9-83a23236046a
begin #from SimLib, changes made
	struct PositionDataDescriptor
	    geometry::Symbol
	    system_size::Int
	    shots::Int
	    ρ::Float64
	    function PositionDataDescriptor(geom, system_size, shots, ρ)
	        geom ∈ GEOMETRIES || error("Unknown geometry: $geom")
	        new(geom, system_size, shots, ρ)
	    end
	end
	
	struct PositionData
	    descriptor::PositionDataDescriptor
	    # coordinates arranged by [N, shot]
	    data::Array{Float64, 2}
		distances::Array{Float64,3} #distance matrix for each shot
	end

	PositionData(desc::PositionDataDescriptor) = PositionData(desc, zeros(Float64, desc.system_size, desc.shots),zeros(Float64,desc.system_size,desc.system_size,desc.shots))
	
	PositionData(args...; kwargs...) = PositionData(PositionDataDescriptor(args..., kwargs...))
end

# ╔═╡ 18643ae4-bd94-451d-a11f-a84b00e7b258
begin #from SimLib, changes made
	function create_positions!(empty_posdata; fail_rate=0.5)
	    logmsg("Generating positions for $(empty_posdata.descriptor)")
	    N = empty_posdata.descriptor.system_size
	    shots = empty_posdata.descriptor.shots
		ρ = empty_posdata.descriptor.ρ
	
	    max_misses = shots*fail_rate
	    sampler = geometry_from_density(empty_posdata.descriptor.geometry, ρ, N, 1) #dim
	    misses = 0
	    for j in 1:shots
	        positions = Vector{Vector{Float64}}()
	        while length(positions) == 0 && misses < max_misses
	            try
	                positions = sample_blockaded(sampler, N)
	            catch e;
	                misses += 1
	                if !(e isa ErrorException)
	                    logmsg("Got unexpected Error: $e\nCurrent ρ=$ρ, shot=$j, N=$N")
	                end
	            end
	        end
	        if misses >= max_misses
	        	error("sampler $sampler did fail to converge too often (fail rate $(round(fail_rate*100; digits=1))%! Current ρ=$ρ, shot=$j, N=$N")
	        end
	        positions = sort!(positions)
	        empty_posdata.data[:,j] = hcat(positions...) #necessary due to :box output
	    	
			empty_posdata.distances[:,:,j] = distance_matrix(sampler,empty_posdata.data[:,j])
		end
	    empty_posdata # now filled
	end
	
	create(desc::PositionDataDescriptor) = create_positions!(PositionData(desc))
end

# ╔═╡ 8021202c-5de9-441f-a9a7-e27004402fee
md"# Test"

# ╔═╡ 0c1b7f13-68ea-4ac9-8e5a-a3ab957862ab
begin
	N = 12
	H = xxz(N,6)
	ψ0 = random_state(N)
	δt = 0.1
	trange = 0:δt:5
end

# ╔═╡ 2f983bed-b3ff-42bd-a7df-b1dde027b08a
desc1 = PositionDataDescriptor(:noisy_chain_pbc,N,10,1)

# ╔═╡ 7b04b3eb-c7f4-4c89-aa39-3d3da0ba12e1
pos1 = create(desc1)

# ╔═╡ 32d4ec59-0879-4750-a53d-d33413ddded2
pos1.data[:,1]

# ╔═╡ 91648919-8ae1-4473-a2ea-e3add45c7ef8
desc2 = PositionDataDescriptor(:box_pbc,N,100,1.2)

# ╔═╡ 680de750-35c6-496e-bf2a-3555e0492859
pos2 = create(desc2)

# ╔═╡ dc25d2fd-5d2f-4019-9647-bd20239cf752
geom = geometry_from_density(:noisy_chain_pbc,1.0,12,1)

# ╔═╡ bf8d87f5-bc44-4763-972b-b087759d2958
md = distance_matrix(geom,pos1.data[:,1])

# ╔═╡ efa25684-bd1d-4cff-b9c4-100e997c1441
md[:,1]

# ╔═╡ 1c817d6c-def0-4d7c-a128-fa03d9d96bba
interaction = PowerLaw(6)

# ╔═╡ 9a18716e-83da-4f19-9f2e-8afd40bf7330
interaction_matrix(interaction,md)

# ╔═╡ 6323ac55-9bba-48b2-8dc4-4862dca83ecb
H2 = xxz(interaction_matrix(interaction,md))

# ╔═╡ 20c2b783-d430-499f-81fb-cff676c07e5e
H2[1:12,1:12]

# ╔═╡ 598f6758-79a1-406b-8cc4-222b59673c53
begin
	A = LightCones.single_spin_op(LightCones.σz,5,N)
	B = LightCones.single_spin_op(LightCones.σz,8,N)
end

# ╔═╡ 93e8222d-2305-4a61-8733-4979e0cfb63f
plot(trange,2*ones(length(trange))-2*otoc(H,A,B,trange,ψ0),ylims=(0,1.5))

# ╔═╡ 62086965-46f6-46de-8a19-b1e1a45c1bda
plot(trange,2*ones(length(trange))-2*otoc(H2,A,B,trange,ψ0),title="ρ=1",ylims=(0,1.5))

# ╔═╡ Cell order:
# ╠═e0bdd126-9e5f-11ec-1174-e7a679707f91
# ╠═452127a7-c67f-4edd-a550-b87ee1392f5b
# ╠═e4052cfb-5950-4b23-a6c8-0e81ee531fa7
# ╠═bf8244ec-4cab-47a0-9481-15974ed7fd96
# ╠═3baa2897-8c0a-48bf-a36d-e99b70723c43
# ╠═74e362cd-fca3-4e54-9e3c-06f3e6f8ea3a
# ╠═effee07f-d0fa-4bee-a741-62ac8339fefe
# ╠═66cccf86-5d37-4e38-93d9-83a23236046a
# ╠═18643ae4-bd94-451d-a11f-a84b00e7b258
# ╠═8021202c-5de9-441f-a9a7-e27004402fee
# ╠═0c1b7f13-68ea-4ac9-8e5a-a3ab957862ab
# ╠═2f983bed-b3ff-42bd-a7df-b1dde027b08a
# ╠═7b04b3eb-c7f4-4c89-aa39-3d3da0ba12e1
# ╠═32d4ec59-0879-4750-a53d-d33413ddded2
# ╠═91648919-8ae1-4473-a2ea-e3add45c7ef8
# ╠═680de750-35c6-496e-bf2a-3555e0492859
# ╠═dc25d2fd-5d2f-4019-9647-bd20239cf752
# ╠═bf8d87f5-bc44-4763-972b-b087759d2958
# ╠═efa25684-bd1d-4cff-b9c4-100e997c1441
# ╠═1c817d6c-def0-4d7c-a128-fa03d9d96bba
# ╠═9a18716e-83da-4f19-9f2e-8afd40bf7330
# ╠═6323ac55-9bba-48b2-8dc4-4862dca83ecb
# ╠═20c2b783-d430-499f-81fb-cff676c07e5e
# ╠═598f6758-79a1-406b-8cc4-222b59673c53
# ╠═93e8222d-2305-4a61-8733-4979e0cfb63f
# ╠═62086965-46f6-46de-8a19-b1e1a45c1bda
