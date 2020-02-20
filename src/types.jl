
# Higher level supertype :
abstract type AbstractBehavior{T} end

abstract type Viscosity{T} <: AbstractBehavior{T} end
abstract type Elasticity{T} <: AbstractBehavior{T} end
abstract type Plasticity{T} <: AbstractBehavior{T} end

Maybe(T)=Union{T,Nothing}

struct Rheology{ T<:AbstractFloat,
        V<:Maybe(Viscosity{T}),
        E<:Maybe(Elasticity{T}),
        P<:Maybe(Plasticity{T}) }
    viscosity::V
    elasticity::E
    plasticity::P
end

Rheology(;viscosity, elasticity, plasticity) = Rheology(viscosity, elasticity, plasticity)
