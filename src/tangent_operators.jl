### TANGENT OPERATORS ###
Dᵉ_func(i,j,k,l,G,λ) = λ*(δ(i,j)*δ(k,l)) + G*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))

get_elastic_stiffness_tensor(e::Elasticity{Float64}) = SymmetricTensor{4, 3}( (i,j,k,l) -> Dᵉ_func(i,j,k,l,e.G,e.λ))
get_elastic_stiffness_tensor(G,λ) = SymmetricTensor{4, 3}( (i,j,k,l) -> Dᵉ_func(i,j,k,l,G,λ))

const DP = DruckerPrager
"""
    compute_stress_tangent(ϵ, r::Rheology{T,::Nothing,::Nothing,::Elasticity,::DruckerPrager}, s::MaterialState ; elastic_only = false) where {T,D,V<:Nothing,E::Elasticity,P::DruckerPrager}

Computes the consistent tangent operator associated with a Drucker-Prager return mapping to the smooth part of the yield function cone based on current strain `ϵ`. Updates the material states and returns the corrected stress as well as the fourth order consistent tangent tensorial operator.
"""
function compute_stress_tangent(ϵ,
                                r::Rheology{T,Nothing,Nothing,E,P},
                                s::MS,
                                clock::Clock;
                                elastic_only = false) where {T,E<:Elasticity,P<:DruckerPrager,MS<:MaterialState}
    ##### unpack some material parameters
    # elastic
    ν = r.elasticity.ν
    G = r.elasticity.G
    K = r.elasticity.K
    Dᵉ = r.elasticity.Dᵉ
    # plastic
    C = r.plasticity.C
    ϕ = r.plasticity.ϕ
    ψ = r.plasticity.ψ
    H = r.plasticity.H
    ηᵛᵖ = r.plasticity.ηᵛᵖ
    η = r.plasticity.η
    η̅ = r.plasticity.η̅
    ξ = r.plasticity.ξ

    ##### Evaluation of trial values
    ϵᵉ_trial = ϵ - s.ϵᵖ # trial-elastic_strain
    σ_trial = Dᵉ ⊡ ϵᵉ_trial # trial-stress
    if elastic_only
        s.temp_ϵ = ϵ
        s.temp_σ = σ_trial
        return s.temp_σ, Dᵉ
    end
    p_trial = 1/3 * tr(σ_trial) # trial pressure, negative in compression
    s_trial = dev(σ_trial) # trial deviatoric stress
    τ_trial = get_τ(s_trial,r.plasticity) # effetive trial-stress (2nd invariant of the deviatoric stress tensor)
    τ_yield = -η*p_trial + ξ*(C + H*s.ϵ̅ᵖ)
    F_trial  = τ_trial - τ_yield # Trial-value of the yield surface

    if F_trial < 0.0 # elastic loading
        s.temp_ϵ = ϵ
        s.temp_σ = σ_trial
        return s.temp_σ, Dᵉ
    else # plastic loading

        # Plastic flow potential gradient
        I2D = SymmetricTensor{2,3}(δ) # second order identity tensor
        ∇Q = ((1/(2*τ_trial)) * s_trial + η̅/3 * I2D)
        ∇Q_scalar = sqrt(2/3 * ∇Q ⊡ ∇Q)

        ##### Compute incremental plastic multiplier
        Δλ_factor = 1/(G + K*η*η̅ + ηᵛᵖ/clock.Δt + ξ*H*∇Q_scalar)
        Δλ = Δλ_factor * F_trial

        ##### Compute unit deviatoric flow vector
        ϵᵉdev_trial = dev(ϵᵉ_trial) #ϵᵉ_trial .- [ϵᵉvol_d3_trial, ϵᵉvol_d3_trial, ϵᵉvol_d3_trial, 0.0]
        norm_ϵᵉdev_t = norm(ϵᵉdev_trial)#sqrt(ϵᵉdev_trial[1]^2 + ϵᵉdev_trial[2]^2 + ϵᵉdev_trial[3]^2 + 2*ϵᵉdev_trial[4]^2) # get the norm of the deviatoric elastic trial strain

        # prevent division by zero if ϵᵉdev_trial is zero valued
        if norm_ϵᵉdev_t != 0.0
            inv_norm = 1.0/norm_ϵᵉdev_t
        else
            inv_norm = 0.0
        end
        uni_dev = ϵᵉdev_trial*inv_norm # unit deviatoric flow vector

        ##### assemble tangent
        Isymdev = SymmetricTensor{4,3}(Isymdev_func) # fourth order deviatoric symmetric tensor
        I2D = SymmetricTensor{2,3}(δ) # second order identity tensor
        A = Δλ_factor
        A_fact = 2G * (1.0 - Δλ/(sqrt(2)*norm_ϵᵉdev_t))
        A_factd3 = A_fact/3
        B_fact = 2G * (Δλ/(sqrt(2)*norm_ϵᵉdev_t) - G*A)
        C_fact = -sqrt(2)*G*A*K
        D_fact = K*(1.0 - K*η*η̅*A)
        D = A_fact * Isymdev +
            B_fact * uni_dev ⊗ uni_dev +
            C_fact * (η * uni_dev⊗I2D + η̅ * I2D⊗uni_dev) +
            D_fact * I2D ⊗ I2D #TODO check whether this form of the D term or the one from the book's code is the good one => (D_fact - A_factd3) or just D_fact.

        Δϵᵖ = Δλ*∇Q # plastic strain increment
        Δϵ̅ᵖ = Δλ*∇Q_scalar # accumulated plastic strain increment
        s.temp_ϵ = ϵ
        s.temp_ϵᵖ = s.ϵᵖ + Δϵᵖ # plastic strain
        s.temp_ϵ̅ᵖ = s.ϵ̅ᵖ + Δϵ̅ᵖ # accumulated plastic strain
        s.temp_σ = σ_trial - Δλ * Dᵉ⊡∇Q
        return s.temp_σ, D
    end
end
