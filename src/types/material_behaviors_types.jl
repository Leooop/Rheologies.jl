### BEHAVIORS ###
# Higher level supertype
# Keyword argument type
"""
`AbstractBehavior` is the abstract supertype of concrete types that contain mechanical properties of the material.
Its subtypes are `Damage`, `Viscosity`, `Elasticity` and `Plasticity`.

These concrete subtypes can be defined with functional fields of a position vector in order to create a rheology that can vary over space, or with a float scalar if the property is constant in space.
This variability or constantness attribute is expressed via the type parameter of `AbstractBehavior` concrete subtypes like `Elasticity{T}` where T can be `Function` or `Float64`. This `Function` type parameter allows the construction of an array containing a `Rheology{Float64}` for each cell with the right scalar properties depending on the cell location.
"""
abstract type AbstractBehavior{T<:F64orFunc} end

# first order abstract hierarchy
abstract type Damage{T} <: AbstractBehavior{T} end
abstract type Plasticity{T} <: AbstractBehavior{T} end

# second order abstract hierarchy
#abstract type ViscoPlasticity{T} <: Plasticity{T} end

### VISCOSITY ###
"""
A `Viscosity` instance contains a viscosity `η`.
"""
@kwdef struct Viscosity{T,T1} <: AbstractBehavior{T}
    η::T1
    function Viscosity(η)
        if η isa Function
            T = Function
        elseif η isa AbstractFloat
            T = typeof(η)
        else
            T = Float64
            η = Float64(η)
        end
        new{T,typeof(η)}(η)
    end
end
Viscosity{T,T1}(η) where {T,T1} = Viscosity(η)

## Elastic concrete types :
# @kwdef struct IncompressibleElasticity{T} <: Elasticity{T}
#     E::Maybe(F64orFunc) = nothing
#     ν::Float64 = 0.4999999
#     G::Maybe(F64orFunc) = nothing
#     K::Maybe(F64orFunc) = nothing
#     λ::Maybe(F64orFunc) = nothing
#     Dᵉ::Maybe(Tensor{4,3}) = nothing
#     function IncompressibleElasticity(E,ν,G,K,λ,Dᵉ)
#         if  (ν >= 0.50) || (ν < 0.49)
#             ν = 0.4999999
#             @warn "Poisson ratio can not be lower than 0.49 and >=0.5 to ensure that the material is close to incompressible. Set to 0.4999999 by default"
#         end
#         T2 = any(isa.((E,G,K,λ),Ref(Function))) ? Function : Float64
#         new{T2}(E,ν,G,K,λ,Dᵉ)
#     end
# end
# IncompressibleElasticity{T}(E,ν,G,K,λ,Dᵉ) where {T} = IncompressibleElasticity(E,ν,G,K,λ,Dᵉ)
#
# function IncompressibleElasticity(E::T,ν::Float64,G::Nothing,K::Nothing,λ::Nothing,Dᵉ::Nothing) where {T<:F64orFunc}
#     if T <: Function
#         return IncompressibleElasticity(E, ν, x -> G_from_Eν(E(x),ν), x -> K_from_Eν(E(x),ν), x -> λ_from_Eν(E(x),ν), nothing)
#     else
#         G = G_from_Eν(E,ν)
#         λ = λ_from_Eν(E,ν)
#         Dᵉ = get_elastic_stiffness_tensor(G,λ)
#         return IncompressibleElasticity(E, ν, G, K_from_Eν(E,ν), λ, Dᵉ)
#     end
# end
#
# function IncompressibleElasticity(E::Nothing,ν::Float64,G::T,K::Nothing,λ::Nothing,Dᵉ::Nothing) where {T<:F64orFunc}
#     if T <: Function
#         return IncompressibleElasticity(x -> E_from_Gν(G(x),ν), ν, G, x -> K_from_Gν(G(x),ν), x -> λ_from_Gν(G(x),ν), nothing)
#     else
#         λ = λ_from_Gν(G,ν)
#         Dᵉ = get_elastic_stiffness_tensor(G,λ)
#         return IncompressibleElasticity(E_from_Gν(G,ν), ν, G, K_from_Gν(G,ν), λ, Dᵉ)
#     end
# end
#
# function IncompressibleElasticity(E::Nothing,ν::Float64,G::Nothing,K::T,λ::Nothing,Dᵉ::Nothing) where {T<:F64orFunc}
#     if T <: Function
#         return IncompressibleElasticity(x -> E_from_Kν(K(x),ν), ν, x -> G_from_Kν(K(x),ν) , K, x -> λ_from_Kν(K(x),ν), nothing)
#     else
#         G = G_from_Kν(K,ν)
#         λ = λ_from_Kν(K,ν)
#         Dᵉ = get_elastic_stiffness_tensor(G,λ)
#         return IncompressibleElasticity(E_from_Kν(K,ν), ν, G, K, λ, Dᵉ)
#     end
# end

"""
An `Elasticity` instance contains elastic properties : moduli `E`, `ν`, `G`, `K`, `λ` and the associated stiffness tensor `Dᵉ`.
"""
@kwdef struct Elasticity{T,T1,T2,T3,T4,T5,T6} <: AbstractBehavior{T}
    E::T1 = nothing
    ν::T2 = nothing
    G::T3 = nothing
    K::T4 = nothing
    λ::T5 = nothing
    Dᵉ::T6 = nothing
    function Elasticity(E,ν,G,K,λ,Dᵉ)
        T = any(isa.((E,ν,G,K,λ),Function)) ? Function : Float64
        new{T,typeof(E),typeof(ν),typeof(G),typeof(K),typeof(λ),typeof(Dᵉ)}(E,ν,G,K,λ,Dᵉ)
    end
end
Elasticity{T,T1,T2,T3,T4,T5,T6}(E,ν,G,K,λ,Dᵉ) where {T,T1,T2,T3,T4,T5,T6} = Elasticity(E,ν,G,K,λ,Dᵉ)

"""
    Elasticity(; <keyword arguments>)
Construct an `Elasticity` instance by giving a pair of keywords arguments, either `E` and `ν` or `G` and `K`.
"""
function Elasticity(E::F64orFunc,ν::F64orFunc,G::Nothing,K::Nothing,λ::Nothing,Dᵉ::Nothing)
    if (E isa Function) & (ν isa Function)
        return Elasticity(E, ν, x -> G_from_Eν(E(x),ν(x)), x -> K_from_Eν(E(x),ν(x)), x -> λ_from_Eν(E(x),ν(x)), nothing)
    elseif E isa Function
        return Elasticity(E, ν, x -> G_from_Eν(E(x),ν), x -> K_from_Eν(E(x),ν), x -> λ_from_Eν(E(x),ν), nothing)
    elseif ν isa Function
        return Elasticity(E, ν, x -> G_from_Eν(E,ν(x)), x -> K_from_Eν(E,ν(x)), x -> λ_from_Eν(E,ν(x)), nothing)
    else
        G = G_from_Eν(E,ν)
        λ = λ_from_Eν(E,ν)
        Dᵉ = get_elastic_stiffness_tensor(G,λ)
        return Elasticity(E, ν, G, K_from_Eν(E,ν), λ, Dᵉ)
    end
end

function Elasticity(E::Nothing,ν::Nothing,G::F64orFunc,K::F64orFunc,λ::Nothing,Dᵉ::Nothing)
    if (G isa Function) & (K isa Function)
        return Elasticity(x -> E_from_KG(K(x),G(x)), x -> ν_from_KG(K(x),G(x)), G, K, x -> λ_from_KG(K(x),G(x)), nothing)
    elseif G isa Function
        return Elasticity(x -> E_from_KG(K,G(x)), x -> ν_from_KG(K,G(x)), G, K, x -> λ_from_KG(K,G(x)), nothing)
    elseif K isa Function
        return Elasticity(x -> E_from_KG(K(x),G), x -> ν_from_KG(K(x),G), G, K, x -> λ_from_KG(K(x),G), nothing)
    else
        λ = λ_from_KG(K,G)
        Dᵉ = get_elastic_stiffness_tensor(G,λ)
        return Elasticity(E_from_KG(K,G), ν_from_KG(K,G), G, K, λ, Dᵉ)
    end
end

function Base.show(io::IO, ::MIME"text/plain", elast::E) where{E<:Elasticity}
    Dᵉ_print = (elast.Dᵉ isa Nothing) ? nothing : "fouth order stiffness tensor"
    print(io, "$(typeof(elast))\n",
    "\t├── E : $(elast.E)\n",
    "\t├── ν : $(elast.ν)\n",
    "\t├── G : $(elast.G)\n",
    "\t├── K : $(elast.K)\n",
    "\t├── λ : $(elast.λ)\n",
    "\t└── Dᵉ : $(Dᵉ_print)\n")
end

"""
`DruckerPrager' type. Its instance contains general plastic properties.

# Fields
μ : friction coefficient
ϕ : friction angle
ψ : dilatancy angle
C : cohesion
H : linear hardening modulus
ηᵛᵖ : Kelvin viscoplastic viscosity
mohr_coulomb_approx::Symbol : Drucker-Prager approximation to Mohr-Coulomb yield criterion
η : pressure prefactor of the yield function
η̅ : pressure prefactor of the plastic flow potential
ξ : cohesion prefactor of the yield function
"""
struct DruckerPrager{T,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10} <: Plasticity{T}
    μ::T1
    ϕ::T2
    ψ::T3
    C::T4
    H::T5
    ηᵛᵖ::T6
    mohr_coulomb_approx::T7
    η::T8
    η̅::T9
    ξ::T10
    function DruckerPrager(μ,ϕ,ψ,C,H,ηᵛᵖ,mohr_coulomb_approx)
        T = any(isa.((μ,ϕ,ψ,C,H,ηᵛᵖ),Function)) ? Function : Float64
        η, η̅, ξ = mohr_coulomb_approx_params(mohr_coulomb_approx, ϕ, ψ)
        return new{T,typeof(μ),typeof(ϕ),typeof(ψ),typeof(C),typeof(H),typeof(ηᵛᵖ),typeof(mohr_coulomb_approx),typeof(η), typeof(η̅), typeof(ξ)}(μ,ϕ,ψ,C,H,ηᵛᵖ,mohr_coulomb_approx,η, η̅, ξ)
    end
end


DruckerPrager{T,T1,T2,T3,T4,T5,T6,T7}(μ,ϕ,ψ,C,H,ηᵛᵖ,mc_approx) where {T,T1,T2,T3,T4,T5,T6,T7} = DruckerPrager(μ,ϕ,ψ,C,H,ηᵛᵖ,mc_approx)
DruckerPrager(μ,ϕ,ψ,C,H,ηᵛᵖ,mohr_coulomb_approx,η, η̅, ξ) = DruckerPrager(μ,ϕ,ψ,C,H,ηᵛᵖ,mohr_coulomb_approx)

"""
    DruckerPrager(; <keyword arguments>)

Construct an `DruckerPrager` instance. At least `μ` or `ϕ` have to be given, other parameters are optional.

# Keyword arguments
μ : friction coefficient
ϕ : friction angle
ψ = ϕ : dilatancy angle
C = 0.0 : cohesion
H = 0.0 : linear hardening modulus
ηᵛᵖ = 0.0 : Kelvin viscoplastic viscosity
mohr_coulomb_approx = :planestrain : Drucker-Prager approximation to Mohr-Coulomb yield criterion
"""
DruckerPrager(; μ=nothing, ϕ=nothing, ψ=nothing,C=0.0,H=0.0,ηᵛᵖ=0.0,mc_approx=:planestrain) = DruckerPrager(μ,ϕ,ψ,C,H,ηᵛᵖ,mc_approx)

# Default case returns an error
DruckerPrager(μ::Nothing,ϕ::Nothing,ψ::Nothing,C,H,ηᵛᵖ,mc_approx) = @error "at least a friction coefficient (μ) or angle (ϕ) must be supplied to create a DruckerPrager instance"

# associated cases
DruckerPrager(μ::Float64,ϕ::Nothing,ψ::Nothing,C,H,ηᵛᵖ,mc_approx) = DruckerPrager(μ,atand(μ),atand(μ),C,H,ηᵛᵖ,mc_approx)
DruckerPrager(μ::Function,ϕ::Nothing,ψ::Nothing,C,H,ηᵛᵖ,mc_approx) = DruckerPrager(μ,x->atand(μ(x)),x->atand(μ(x)),C,H,ηᵛᵖ,mc_approx)
DruckerPrager(μ::Nothing,ϕ::Float64,ψ::Nothing,C,H,ηᵛᵖ,mc_approx) = DruckerPrager(tand(ϕ),ϕ,ϕ,C,H,ηᵛᵖ,mc_approx)
DruckerPrager(μ::Nothing,ϕ::Function,ψ::Nothing,C,H,ηᵛᵖ,mc_approx) = DruckerPrager(x->tand(ϕ(x)),ϕ,ϕ,C,H,ηᵛᵖ,mc_approx)

# with dilatancy angle specified
DruckerPrager(μ::Float64,ϕ::Nothing,ψ,C,H,ηᵛᵖ,mc_approx) = DruckerPrager(μ,atand(μ),ψ,C,H,ηᵛᵖ,mc_approx)
DruckerPrager(μ::Function,ϕ::Nothing,ψ,C,H,ηᵛᵖ,mc_approx) = DruckerPrager(μ,x->atand(μ(x)),ψ,C,H,ηᵛᵖ,mc_approx)
DruckerPrager(μ::Nothing,ϕ::Float64,ψ,C,H,ηᵛᵖ,mc_approx) = DruckerPrager(tand(ϕ),ϕ,ψ,C,H,ηᵛᵖ,mc_approx)
DruckerPrager(μ::Nothing,ϕ::Function,ψ,C,H,ηᵛᵖ,mc_approx) = DruckerPrager(x->tand(ϕ(x)),ϕ,ψ,C,H,ηᵛᵖ,mc_approx)

function Base.show(io::IO, ::MIME"text/plain",
                   plast::DruckerPrager)
    print(io, "$(typeof(plast))\n",
    "\t├── μ (friction coefficient)                 : $(plast.μ)\n",
    "\t├── ϕ (friction angle)                       : $(plast.ϕ)\n",
    "\t├── ψ (dilatancy angle)                      : $(plast.ψ)\n",
    "\t├── C (cohesion)                             : $(plast.C)\n",
    "\t├── H (hardening modulus)                    : $(plast.H)\n",
    "\t├── ηᵛᵖ (viscosity)                          : $(plast.ηᵛᵖ)\n",
    "\t├── mohr_coulomb_approx                      : $(plast.mohr_coulomb_approx)\n",
    "\t├── η (yield pressure prefactor)             : $(plast.η)\n",
    "\t├── η̅ (plastic potential pressure prefactor) : $(plast.η̅)\n",
    "\t└── ξ (cohesion prefactor)                   : $(plast.ξ)\n")
end



# @kwdef mutable struct ViscoDruckerPrager{T} <: ViscoPlasticity{T}
#     μ::F64orFunc = 0.0
#     ϕ::F64orFunc = 0.0
#     ψ::F64orFunc = 0.0
#     C::F64orFunc = 0.0
#     H::F64orFunc = 0.0
#     ηᵛᵖ::F64orFunc = 0.0
#     mohr_coulomb_approx::Symbol = :planestrain
#     η::F64orFunc = nothing
#     η̅::F64orFunc = nothing
#     ξ::F64orFunc = nothing
#     function ViscoDruckerPrager(μ,ϕ,ψ,C,H,ηᵛᵖ,mohr_coulomb_approx)
#         T2 = any(isa.((μ,ϕ,ψ,C,H,ηᵛᵖ),Function)) ? Function : typeof(μ)
#         if ϕ == 0
#             if μ != 0
#                 if μ isa Function
#                     return new{T2}(μ, x->atand(μ(x)), ψ, C, H, ηᵛᵖ)
#                 elseif μ isa Float64
#                     return new{T2}(μ, atand(μ), ψ, C, H, ηᵛᵖ)
#                 else
#                     @error "Unexpected type of field μ"
#                 end
#             else
#                 @error "At least ϕ or μ must be provided"
#             end
#         else
#             #μ != 0 && @error "Ambiguity between μ and ϕ, please provide only one of them"
#             if ϕ isa Function
#                 return new{T2}(x->tand(ϕ(x)), ϕ, ψ, C, H, ηᵛᵖ)
#             elseif ϕ isa Float64
#                 return new{T2}(tand(ϕ), ϕ, ψ, C, H, ηᵛᵖ)
#             else
#                 @error "Unexpected type of field ϕ"
#             end
#         end
#     end
# end
# ViscoDruckerPrager{T}(μ,ϕ,ψ,C,H,ηᵛᵖ) where {T} = ViscoDruckerPrager(μ,ϕ,ψ,C,H,ηᵛᵖ)


@kwdef mutable struct VonMises{T,T1} <: Plasticity{T}
    C::T1
    function VonMises(C)
        T = isa(C,Function) ? Function : typeof(C)
        new{T,typeof(C)}(C)
    end
end
VonMises{T,T1}(C) where {T,T1} = VonMises(C)

function Base.show(io::IO, ::MIME"text/plain",
                   plast::VonMises)
    print(io, "$(typeof(plast))\n",
    "\t└── C : $(plast.C)")
end

## Damage concrete types :
### BHAT ROSAKIS SAMMIS micromechanical parameters ###
@kwdef mutable struct BRSDamage{T} <: Damage{T}
    μ::F64orFunc = 0.6 # Friction coef
    β::F64orFunc = 0.1 # Correction factor
    K₁c::F64orFunc = 1.74e6 # Critical stress intensity factor (Pa.m^(1/2))
    a::F64orFunc # Initial flaw size (m)
    ψ::F64orFunc # crack angle to the principal stress (radians)
    D₀::F64orFunc # Initial flaw density
    n::F64orFunc = 34.0 # Stress corrosion index
    l̇₀::F64orFunc = 0.24 # Ref. crack growth rate (m/s)
    H::F64orFunc = 50e3 # Activation enthalpy (J/mol)
    A::F64orFunc = 5.71 # Preexponential factor (m/s)
    function BRSDamage(μ,β,K₁c,a,ψ,D₀,n,l̇₀,H,A)
        T2 = any(isa.((μ,β,K₁c,a,ψ,D₀,n,l̇₀,H,A),Function)) ? Function : typeof(β)
        new{T2}(μ,β,K₁c,a,ψ,D₀,n,l̇₀,H,A)
    end
end
BRSDamage{T}(μ,β,K₁c,a,ψ,D₀,n,l̇₀,H,A) where {T} = BRSDamage(μ,β,K₁c,a,ψ,D₀,n,l̇₀,H,A)



## UPDATE PARAMETERS :
function updatetypeparameter(::Viscosity,p)
    return Viscosity(p...)
end
function updatetypeparameter(::Elasticity,p)
    return Elasticity(p...)
end
function updatetypeparameter(::DruckerPrager,p)
    return DruckerPrager(p...)
end
function updatetypeparameter(::VonMises,p)
    return VonMises(p...)
end
function updatetypeparameter(::BRSDamage,p)
    return BRSDamage(p...)
end
# function updatetypeparameter(::ViscoDruckerPrager,p)
#     return ViscoDruckerPrager(p...)
# end

# would be gould to generate these functions for each Behavior <: AbstractBehavior using code generation, but the following raises an UndefVarError "CompressibleElasticity not defined"

# for behavior in subtypes(AbstractBehavior)
#     for concrete_behavior in subtypes(behavior)
#         typ = Symbol(concrete_behavior)
#         println(typ)
#         eval(
#             :(function updatetypeparameter(b::$(typ))
#                 properties = getproperty.(Ref(b),propertynames(b))
#                 return $(typ)(properties...)
#             end)
#         )
#     end
# end

function mohr_coulomb_approx_params(mohr_coulomb_approx, ϕ, ψ)
    if mohr_coulomb_approx == :planestrain
        if ϕ isa Function
            η = x -> 3*tand(ϕ(x))/(sqrt(9 + 12*tand(ϕ(x))^2))
            ξ = x -> 3/(sqrt(9 + 12*tand(ϕ(x))^2))
        else
            η = 3*tand(ϕ)/(sqrt(9 + 12*tand(ϕ)^2))
            ξ = 3/(sqrt(9 + 12*tand(ϕ)^2))
        end
        if ψ isa Function
            η̅ = x -> 3*tand(ψ(x))/(sqrt(9 + 12*tand(ψ(x))^2))
        else
            η̅ = 3*tand(ψ)/(sqrt(9 + 12*tand(ψ)^2))
        end
    elseif mohr_coulomb_approx == :outeredges
        if ϕ isa Function
            η = x -> 6*sind(ϕ(x))/(sqrt(3)*(3 - sind(ϕ(x))))
            ξ = x -> 6*cosd(ϕ(x))/(sqrt(3)*(3 - sind(ϕ(x))))
        else
            η = 6*sind(ϕ)/(sqrt(3)*(3 - sind(ϕ)))
            ξ = 6*cosd(ϕ)/(sqrt(3)*(3 - sind(ϕ)))
        end
        if ψ isa Function
            η̅ = x -> 6*sind(ψ(x))/(sqrt(3)*(3 - sind(ψ(x))))
        else
            η̅ = 6*sind(ψ)/(sqrt(3)*(3 - sind(ψ)))
        end
    elseif mohr_coulomb_approx == :inneredges
        if ϕ isa Function
            η = x -> 6*sind(ϕ(x))/(sqrt(3)*(3 + sind(ϕ(x))))
            ξ = x -> 6*cosd(ϕ(x))/(sqrt(3)*(3 + sind(ϕ(x))))
        else
            η = 6*sind(ϕ)/(sqrt(3)*(3 + sind(ϕ)))
            ξ = 6*cosd(ϕ)/(sqrt(3)*(3 + sind(ϕ)))
        end
        if ψ isa Function
            η̅ = x -> 6*sind(ψ(x))/(sqrt(3)*(3 + sind(ψ(x))))
        else
            η̅ = 6*sind(ψ)/(sqrt(3)*(3 + sind(ψ)))
        end
    else
        @error "Drucker Prager mohr coulomb approximation $(mohr_coulomb_approx) is not available. Use :planestrain, :outeredges or :inner edges instead"
    end
    return η, η̅, ξ
end
