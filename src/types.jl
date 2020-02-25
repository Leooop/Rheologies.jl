
using Base: @kwdef

### BEHAVIORS ###
# Higher level supertype :
abstract type AbstractBehavior{T} end

abstract type Viscosity{T} <: AbstractBehavior{T} end
abstract type Elasticity{T} <: AbstractBehavior{T} end
abstract type Plasticity{T} <: AbstractBehavior{T} end

@kwdef struct IncompressibleElasticity{T} <: Elasticity{T}
    E::T
    ν::T
    function IncompressibleElasticity(E,ν)
        if ν <= 0.49
            ν = 0.4999999
            @warn "Poisson ratio can not be lower than 0.49 to ensure that the material is
            close to incompressible. Set to 0.4999999 by default"
        end
        new{typeof(E)}(E,convert(typeof(E),ν))
    end
end

@kwdef struct CompressibleElasticity{T} <: Elasticity{T}
    E::T
    ν::T
end

# function Base.show(io::IO, ::MIME"text/plain", elast::E) where{E<:Elasticity}
#     println(io, "Elasticity of type $(typeof(elast))")
#     println("    -> E : $(elast.E)")
#     println("    -> ν : $(elast.ν)")
# end

get_G(e::Elasticity) = e.E/(2*(1+e.ν))
get_K(e::Elasticity) = e.E/(3*(1-2*e.ν))
## RHEOLOGY ###
Maybe(T)=Union{T,Nothing}

struct Rheology{ T<:AbstractFloat,
        V<:Maybe(Viscosity{T}),
        E<:Maybe(Elasticity{T}),
        P<:Maybe(Plasticity{T}) }
    viscosity::V
    elasticity::E
    plasticity::P
end

Rheology(T ;viscosity, elasticity, plasticity) = Rheology(T,viscosity, elasticity, plasticity)

function Base.show(io::IO, ::MIME"text/plain", rheology::Rheology)
    println(io, "Rheology type")
    println("    -> viscosity  : $(typeof(rheology.viscosity))")
    println("    -> elasticity : $(typeof(rheology.elasticity))")
    println("    -> plasticity : $(typeof(rheology.plasticity))")
end

rheology_summary(r::Rheology{T,Nothing,E,Nothing}) where{T,E<:Elasticity} = "elastic"
rheology_summary(r::Rheology{T,V,E,Nothing}) where{T,V<:Viscosity,E<:Elasticity} = "visco-elastic"
rheology_summary(r::Rheology{T,V,Nothing,Nothing}) where{T,V<:Viscosity} = "viscous"
rheology_summary(r::Rheology{T,Nothing,E,P}) where{T,E<:Elasticity,P<:Plasticity} = "elasto-plastic"
rheology_summary(r::Rheology{T,V,E,P}) where{T,V<:Viscosity,E<:Elasticity,P<:Plasticity} = "visco-elasto-plastic"
### BOUNDARY CONDITIONS ###

@kwdef struct BoundaryConditions
    dirichlet::Dict
    neumann::Maybe(Dict)
end
