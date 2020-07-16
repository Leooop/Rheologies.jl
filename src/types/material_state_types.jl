### MATERIAL STATES ###
"""
`MaterialState` is the abstract supertype of objets storing the stresses and strains states of the material at the location of one quadrature point.
The states of the whole material are then stored in a `Vector{Vector{<:MaterialState}} (`MaterialState`s contained in a cell, and the cell in a grid).
Stores stress and strain informations as tensorial quantities.
"""
abstract type MaterialState{S <: SecondOrderTensor} end


function update_material_state!(model)
    for cell_states in model.material_state
        foreach(update_state!, cell_states)
    end
end

function copy_temp_state!(model,mp)
    for (cell_id, cell_states) in enumerate(model.material_state)
        for (qp_id,state) in enumerate(cell_states)
            for field in propertynames(state)
                str_field = String(field)
                if length(str_field) > 4
                    (str_field[1:4] == "temp") && setproperty!(state,field,getproperty(mp[cell_id][qp_id],field))
                end
            end
        end
    end
end
##### Elastic or viscous #####
"""
`BasicMaterialState` is used for visco-elastic materials. It stores the elastic strains and stresses.

# Fields
ϵ : strain tensor
σ : stress tensor
temp_ϵ : temporary strain tensor
temp_σ : temporary stress tensor
"""
mutable struct BasicMaterialState{S} <: MaterialState{S}
    # Store converged values
    ϵ::S # total strain
    σ::S # stress

    # Store temporary values used during equilibrium iterations
    temp_ϵ::S
    temp_σ::S
end

function BasicMaterialState()
    return BasicMaterialState(
                zero(SymmetricTensor{2, 3}),
                zero(SymmetricTensor{2, 3}),
                zero(SymmetricTensor{2, 3}),
                zero(SymmetricTensor{2, 3}))
end

function update_state!(state::BasicMaterialState)
    state.ϵ = state.temp_ϵ
    state.σ = state.temp_σ
end


##### Plastic #####
"""
`PlasticMaterialState` is used for visco-elasto-plastic materials. It stores the elastic strains and stresses as well as the plastic strains and the accumulated plastic strain.

# Fields
ϵᵖ : plastic strain tensor
ϵ̅ᵖ : accumulated plastic strain
ϵ  : strain tensor
σ  : stress tensor
temp_ϵᵖ : temporary plastic strain tensor
temp_ϵ̅ᵖ : temporary accumulated plastic strain
temp_ϵ : temporary strain tensor
temp_σ : temporary stress tensor
"""
mutable struct PlasticMaterialState{S} <: MaterialState{S}
    # Store converged values
    ϵᵖ::S # plastic strain
    ϵ̅ᵖ::Float64 # accumulated plastic strain
    ϵ::S # total strain
    σ::S # stress

    # Store temporary values used during equilibrium iterations
    temp_ϵᵖ::S
    temp_ϵ̅ᵖ::Float64
    temp_ϵ::S
    temp_σ::S
end
function PlasticMaterialState()
    return PlasticMaterialState(
                zero(SymmetricTensor{2, 3}),
                0.0,
                zero(SymmetricTensor{2, 3}),
                zero(SymmetricTensor{2, 3}),
                zero(SymmetricTensor{2, 3}),
                0.0,
                zero(SymmetricTensor{2, 3}),
                zero(SymmetricTensor{2, 3}))
end

function update_state!(state::PlasticMaterialState)
    state.ϵᵖ = state.temp_ϵᵖ
    state.ϵ̅ᵖ = state.temp_ϵ̅ᵖ
    state.ϵ = state.temp_ϵ
    state.σ = state.temp_σ
end

function set_temp_state!(r::Rheology{TT,TD,Nothing,TE,TP},clock,state,σ,ϵ) where {TT,TD,TE,TP}
    Dᵉ = r.elasticity.Dᵉ
    Δϵ = ϵ - state.ϵ
    invDᵉ = get_elastic_compliance_tensor(r)
    Δϵᵖ = Δϵ - invDᵉ⊡(σ - state.σ)
    ϵ̇ᵖ = Δϵᵖ/clock.Δt
    Δϵ̅ᵖ = sqrt((2/3) * ϵ̇ᵖ ⊡ ϵ̇ᵖ) * clock.Δt

    state.temp_ϵᵖ = state.ϵᵖ + Δϵᵖ
    state.temp_ϵ̅ᵖ = state.ϵ̅ᵖ + Δϵ̅ᵖ
    state.temp_σ = σ
    state.temp_ϵ = ϵ
end
##### Damaged plastic #####
"""
`DamagedPlasticMaterialState` is used for damaged visco-elasto-plastic materials.
It stores the elastic strains and stresses as well as the plastic strains and the accumulated plastic strain and a damage state variable.

# Fields
D  : damage state variable
ϵᵖ : plastic strain tensor
ϵ̅ᵖ : accumulated plastic strain
ϵ  : strain tensor
σ  : stress tensor
temp_D  : temporary damage state variable
temp_ϵᵖ : temporary plastic strain tensor
temp_ϵ̅ᵖ : temporary accumulated plastic strain
temp_ϵ  : temporary strain tensor
temp_σ  : temporary stress tensor
"""
mutable struct DamagedPlasticMaterialState{S} <: MaterialState{S}
    # Store converged values
    D::Float64
    ϵᵖ::S # plastic strain
    ϵ̅ᵖ::Float64 # cumulated plastic strain
    ϵ::S # total strain
    σ::S # stress

    # Store temporary values used during equilibrium iterations
    temp_D::Float64
    temp_ϵᵖ::S
    temp_ϵ̅ᵖ::Float64
    temp_ϵ::S
    temp_σ::S
end

function DamagedPlasticMaterialState(r::Rheology)
    return DamagedPlasticMaterialState(
                r.damage.D₀,
                zero(SymmetricTensor{2, 3}),
                0.0,
                zero(SymmetricTensor{2, 3}),
                zero(SymmetricTensor{2, 3}),
                r.damage.D₀,
                zero(SymmetricTensor{2, 3}),
                0.0,
                zero(SymmetricTensor{2, 3}),
                zero(SymmetricTensor{2, 3}))
end


function update_state!(state::DamagedPlasticMaterialState)
    state.D = state.temp_D
    state.ϵᵖ = state.temp_ϵᵖ
    state.ϵ̅ᵖ = state.temp_ϵ̅ᵖ
    state.ϵ = state.temp_ϵ
    state.σ = state.temp_σ
end
