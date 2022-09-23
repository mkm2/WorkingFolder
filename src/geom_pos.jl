module GeomPos

using LinearAlgebra
using SpinSymmetry
using ..LightCones

export GEOMETRIES, VOLUME_COEFFS, SAFE_RHO_RANGES
export PositionData, PositionDataDescriptor
export geometry_from_density, create 


 #from SimLib, no changes made
const GEOMETRIES = [:box, :box_pbc, :noisy_chain, :noisy_chain_pbc, :regular_chain, :regular_chain_pbc]

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


#from SimLib, changes made
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

PositionData(desc::PositionDataDescriptor) = PositionData(desc, zeros(Float64, desc.system_size, desc.shots), zeros(Float64,desc.system_size,desc.system_size,desc.shots))

PositionData(args...; kwargs...) = PositionData(PositionDataDescriptor(args..., kwargs...))

#from SimLib, changes made
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


end #module