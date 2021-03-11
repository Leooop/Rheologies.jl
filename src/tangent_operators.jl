### TANGENT OPERATORS ###

"""
    compute_stress_tangent(ϵ, r::Rheology{T,::Nothing,::Nothing,::Elasticity,::DruckerPrager}, s::MaterialState ; elastic_only = false) where {T,D,V<:Nothing,E::Elasticity,P::DruckerPrager}

Computes the consistent tangent operator associated with a Drucker-Prager return mapping to the smooth part of the yield function cone based on current strain `ϵ`. Updates the material states and returns the corrected stress as well as the fourth order consistent tangent tensorial operator.
"""
function compute_stress_tangent_old(ϵ,
                                r::Rheology{T,Nothing,Nothing,E,P},
                                s::MaterialState,
                                clock::Clock;
                                noplast = false) where {T,E<:Elasticity,P<:DruckerPrager}

    # unpack elastic stiffness tensor
    Dᵉ = r.elasticity.Dᵉ

    ##### Evaluation of trial values
    ϵᵉ_trial = ϵ - s.ϵᵖ # trial-elastic_strain
    σ_trial = Dᵉ ⊡ ϵᵉ_trial # trial-stress

    if noplast # during first newton iteration to homogenize strain field
        s.temp_ϵ = ϵ
        s.temp_σ = σ_trial
        return s.temp_σ, Dᵉ
    end

    # unpack plastic properties
    C = r.plasticity.C
    H = r.plasticity.H
    η = r.plasticity.η
    ξ = r.plasticity.ξ

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
        # unpack
        G = r.elasticity.G
        K = r.elasticity.K
        ηᵛᵖ = r.plasticity.ηᵛᵖ
        η̅ = r.plasticity.η̅

        # Plastic flow potential gradient
        I2D = SymmetricTensor{2,3}(δ) # second order identity tensor
        ∇Q = ((1/(2*τ_trial)) * s_trial + η̅/3 * I2D)
        #∇Q_scalar = sqrt(2/3 * ∇Q ⊡ ∇Q)

        ##### Compute incremental plastic multiplier
        Δλ_factor = 1/(G + K*η*η̅ + ηᵛᵖ/clock.Δt + ξ^2*H) #ξ*H*∇Q_scalar
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
        Δϵ̅ᵖ = Δλ*ξ#Δλ*∇Q_scalar # accumulated plastic strain increment
        s.temp_ϵ = ϵ
        s.temp_ϵᵖ = s.ϵᵖ + Δϵᵖ # plastic strain
        s.temp_ϵ̅ᵖ = s.ϵ̅ᵖ + Δϵ̅ᵖ # accumulated plastic strain
        s.temp_σ = σ_trial - Δλ * Dᵉ⊡∇Q
        return s.temp_σ, D
    end
end

#TODO put back this function if the version iterating over Δλ doesn't work.
function compute_stress_tangent_prev(ϵ,
                                r::Rheology{T,Nothing,Nothing,E,P},
                                s::MaterialState,
                                clock::Clock;
                                noplast = false) where {T,E<:Elasticity,P<:DruckerPrager}

    # unpack elastic stiffness tensor
    Dᵉ = r.elasticity.Dᵉ

    ##### Evaluation of trial values
    ϵᵉ_trial = ϵ - s.ϵᵖ # trial-elastic_strain
    σ_trial = Dᵉ ⊡ ϵᵉ_trial # trial-stress

    if noplast # during first newton iteration to homogenize strain field
        s.temp_ϵ = ϵ
        s.temp_σ = σ_trial
        return s.temp_σ, Dᵉ
    end

    # unpack plastic properties
    C = r.plasticity.C
    H = r.plasticity.H
    η = r.plasticity.η
    ξ = r.plasticity.ξ

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
        # unpack
        G = r.elasticity.G
        K = r.elasticity.K
        ηᵛᵖ = r.plasticity.ηᵛᵖ
        η̅ = r.plasticity.η̅

        # Plastic flow potential gradient
        I2D = SymmetricTensor{2,3}(δ) # second order identity tensor
        ∇Q = ((1/(2*τ_trial)) * s_trial + η̅/3 * I2D)
        #∇Q_scalar = sqrt(2/3 * ∇Q ⊡ ∇Q)

        ##### Compute incremental plastic multiplier
        Δλ_factor = 1/(G + K*η*η̅ + ηᵛᵖ/clock.Δt + ξ^2*H) # + ηᵛᵖ/clock.Δt pour visco-plasticité
        Δλ = Δλ_factor * F_trial

        ##### Test the validity of the return mapping to the smooth portion of the cone
        is_valid = (τ_trial - G*Δλ >= 0)

        ##### Use relevant return mapping algorithm

        if is_valid # return to the smooth portion of the cone

            ##### Compute unit deviatoric flow vector
            ϵᵉdev_trial = dev(ϵᵉ_trial) #ϵᵉ_trial .- [ϵᵉvol_d3_trial, ϵᵉvol_d3_trial, ϵᵉvol_d3_trial, 0.0]
            norm_ϵᵉdev_t = norm(ϵᵉdev_trial)#sqrt(ϵᵉdev_trial[1]^2 + ϵᵉdev_trial[2]^2 + ϵᵉdev_trial[3]^2 + 2*ϵᵉdev_trial[4]^2) # get the norm of the deviatoric elastic trial strain

            # prevent division by zero if ϵᵉdev_trial is zero valued
            if norm_ϵᵉdev_t != 0.0
                inv_norm = 1.0/norm_ϵᵉdev_t
            else
                inv_norm = 0.0
            end
            uni_dev = ϵᵉdev_trial*inv_norm # unit deviatoric flow tensor

            ##### assemble tangent
            Isymdev = SymmetricTensor{4,3}(Isymdev_func) # fourth order deviatoric symmetric tensor
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
            Δϵ̅ᵖ = Δλ*ξ#Δλ*∇Q_scalar # accumulated plastic strain increment
            s.temp_σ = σ_trial - Δλ * Dᵉ⊡∇Q

        else # return to the apex
            α = ξ/η
            β = ξ/η̅
            Δϵᵖ_vol = (p_trial - (C + H*s.ϵ̅ᵖ)*β) / (H*α*β + K)
            Δϵᵖ = 1/3 * Δϵᵖ_vol * I2D
            Δϵ̅ᵖ = α * Δϵᵖ_vol
            s.temp_σ = (p_trial - K*Δϵᵖ_vol)*I2D

            D = K * (1 - K/(K + α*β*H)) * I2D ⊗ I2D
        end

        ### update strains
        s.temp_ϵ = ϵ
        s.temp_ϵᵖ = s.ϵᵖ + Δϵᵖ # plastic strain
        s.temp_ϵ̅ᵖ = s.ϵ̅ᵖ + Δϵ̅ᵖ # accumulated plastic strain

        return s.temp_σ, D
    end
end

function compute_stress_tangent(ϵ,
                                r::Rheology{T,Nothing,Nothing,E,P},
                                s::MaterialState,
                                clock::Clock;
                                noplast = false) where {T,E<:Elasticity,P<:DruckerPrager}

    # unpack elastic stiffness tensor
    Dᵉ = r.elasticity.Dᵉ

    ##### Evaluation of trial values
    ϵᵉ_trial = ϵ - s.ϵᵖ # trial-elastic_strain
    σ_trial = Dᵉ ⊡ ϵᵉ_trial # trial-stress

    if noplast # during first newton iteration to homogenize strain field
        s.temp_ϵ = ϵ
        s.temp_σ = σ_trial
        return s.temp_σ, Dᵉ
    end

    # unpack plastic properties
    C = r.plasticity.C
    H = r.plasticity.H
    η = r.plasticity.η
    ξ = r.plasticity.ξ
    ηᵛᵖ = r.plasticity.ηᵛᵖ
    η̅ = r.plasticity.η̅
    # unpack elastic properties
    G = r.elasticity.G
    K = r.elasticity.K

    p_trial = 1/3 * tr(σ_trial) # trial pressure, negative in compression
    s_trial = dev(σ_trial) # trial deviatoric stress
    τ_trial = get_τ(s_trial,r.plasticity) # effetive trial-stress (2nd invariant of the deviatoric stress tensor)
    τ_yield = -η*p_trial + ξ*(C + H*s.ϵ̅ᵖ)
    F_trial_0  = τ_trial - τ_yield # Trial-value of the yield surface

    converged_Δλ = false
    Δλ₀ = 0.0
    Δλ = 0.0
    Δλ_factor = 0.0
    F_trial = 0.0

    while !converged_Δλ # iterate until Δλ converge
        # compute F_trial with current Δλ₀
        F_trial = F_trial_0 - ηᵛᵖ*Δλ₀/clock.Δt
        ##### Compute incremental plastic multiplier
        Δλ_factor = 1/(G + K*η*η̅ + ηᵛᵖ/clock.Δt + ξ^2*H) # + ηᵛᵖ/clock.Δt pour visco-plasticité
        Δλ = Δλ_factor * F_trial
        # Check error
        (abs(Δλ-Δλ₀)/Δλ₀ <= 1e-8) && (converged_Δλ = true)
        Δλ₀ = Δλ
    end

    if F_trial < 0.0 # elastic loading
        s.temp_ϵ = ϵ
        s.temp_σ = σ_trial
        return s.temp_σ, Dᵉ
    else # plastic loading

        # Plastic flow potential gradient
        I2D = SymmetricTensor{2,3}(δ) # second order identity tensor
        ∇Q = ((1/(2*τ_trial)) * s_trial + η̅/3 * I2D)
        #∇Q_scalar = sqrt(2/3 * ∇Q ⊡ ∇Q)


        ##### Test the validity of the return mapping to the smooth portion of the cone
        is_valid = (τ_trial - G*Δλ >= 0)

        ##### Use relevant return mapping algorithm

        if is_valid # return to the smooth portion of the cone

            ##### Compute unit deviatoric flow vector
            ϵᵉdev_trial = dev(ϵᵉ_trial) #ϵᵉ_trial .- [ϵᵉvol_d3_trial, ϵᵉvol_d3_trial, ϵᵉvol_d3_trial, 0.0]
            norm_ϵᵉdev_t = norm(ϵᵉdev_trial)#sqrt(ϵᵉdev_trial[1]^2 + ϵᵉdev_trial[2]^2 + ϵᵉdev_trial[3]^2 + 2*ϵᵉdev_trial[4]^2) # get the norm of the deviatoric elastic trial strain

            # prevent division by zero if ϵᵉdev_trial is zero valued
            if norm_ϵᵉdev_t != 0.0
                inv_norm = 1.0/norm_ϵᵉdev_t
            else
                inv_norm = 0.0
            end
            uni_dev = ϵᵉdev_trial*inv_norm # unit deviatoric flow tensor

            ##### assemble tangent
            Isymdev = SymmetricTensor{4,3}(Isymdev_func) # fourth order deviatoric symmetric tensor
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
            Δϵ̅ᵖ = Δλ*ξ#Δλ*∇Q_scalar # accumulated plastic strain increment
            s.temp_σ = σ_trial - Δλ * Dᵉ⊡∇Q

        else # return to the apex
            α = ξ/η
            β = ξ/η̅
            Δϵᵖ_vol = (p_trial - (C + H*s.ϵ̅ᵖ)*β) / (H*α*β + K)
            Δϵᵖ = 1/3 * Δϵᵖ_vol * I2D
            Δϵ̅ᵖ = α * Δϵᵖ_vol
            s.temp_σ = (p_trial - K*Δϵᵖ_vol)*I2D

            D = K * (1 - K/(K + α*β*H)) * I2D ⊗ I2D
        end

        ### update strains
        s.temp_ϵ = ϵ
        s.temp_ϵᵖ = s.ϵᵖ + Δϵᵖ # plastic strain
        s.temp_ϵ̅ᵖ = s.ϵ̅ᵖ + Δϵ̅ᵖ # accumulated plastic strain

        return s.temp_σ, D
    end
end

function compute_stress_tangent(ϵ,
                                r::Rheology{T,Nothing,V,E,Nothing},
                                s::MaterialState,
                                clock::Clock;
                                noplast = false) where {T,V<:Viscosity,E<:Elasticity}

    #@info "good stress tangent function"

    ##### viscoelastic tangent operator
    Gᵛᵉ = 1/(1/r.elasticity.G + clock.Δt/r.viscosity.η)
    λᵛᵉ = λ_from_KG(r.elasticity.K,Gᵛᵉ)
    Dᵛᵉ = get_elastic_stiffness_tensor(Gᵛᵉ,λᵛᵉ)

    ##### stress and strain update
    Gratio = Gᵛᵉ/r.elasticity.G
    I2D = SymmetricTensor{2,3}(δ) # second order identity tensor
    σ⁰ = (1/3)*tr(s.σ)*I2D + Gratio*dev(s.σ)
    s.temp_σ = σ⁰ + Dᵛᵉ ⊡ (ϵ - s.ϵ)
    s.temp_ϵ = ϵ

    return s.temp_σ, Dᵛᵉ
end

function compute_stress_tangent(ϵ,
                                r::Rheology{T,Nothing,V,E,P},
                                s::MaterialState,
                                clock::Clock;
                                noplast = false) where {T,V<:Viscosity,E<:Elasticity,P<:DruckerPrager}

    ##### viscoelastic tangent operator
    Gᵛᵉ = 1/(1/r.elasticity.G + clock.Δt/r.viscosity.η)
    λᵛᵉ = λ_from_KG(r.elasticity.K,Gᵛᵉ)
    Dᵛᵉ = get_elastic_stiffness_tensor(Gᵛᵉ,λᵛᵉ)

    # trial stress
    Gratio = Gᵛᵉ/r.elasticity.G
    I2D = SymmetricTensor{2,3}(δ) # second order identity tensor
    σ⁰ = (1/3)*tr(s.σ)*I2D + Gratio*dev(s.σ)
    # TODO find good expression for ϵᵉ_trial
    #ϵᵉ_trial = ϵ - s.ϵᵖ # should be bad because doesn't take viscous deformation into account
    ϵᵉ_trial = ϵ - s.ϵ

    σ_trial = σ⁰ + Dᵛᵉ ⊡ (ϵ - s.ϵ)

    if noplast # during first newton iteration to homogenize strain field
        s.temp_ϵ = ϵ
        s.temp_σ = σ_trial
        return s.temp_σ, Dᵛᵉ
    end

    # unpack plastic properties
    C = r.plasticity.C
    H = r.plasticity.H
    η = r.plasticity.η
    ξ = r.plasticity.ξ

    p_trial = 1/3 * tr(σ_trial) # trial pressure, negative in compression
    s_trial = dev(σ_trial) # trial deviatoric stress
    τ_trial = get_τ(s_trial,r.plasticity) # effetive trial-stress (2nd invariant of the deviatoric stress tensor)
    τ_yield = -η*p_trial + ξ*(C + H*s.ϵ̅ᵖ)
    F_trial  = τ_trial - τ_yield # Trial-value of the yield surface

    if F_trial < 0.0 # elastic loading
        s.temp_ϵ = ϵ
        s.temp_σ = σ_trial
        return s.temp_σ, Dᵛᵉ
    else # plastic loading
        # unpack
        K = r.elasticity.K
        ηᵛᵖ = r.plasticity.ηᵛᵖ
        η̅ = r.plasticity.η̅

        # Plastic flow potential gradient
        ∇Q = ((1/(2*τ_trial)) * s_trial + η̅/3 * I2D)
        #∇Q_scalar = sqrt(2/3 * ∇Q ⊡ ∇Q)

        ##### Compute incremental plastic multiplier
        Δλ_factor = 1/(Gᵛᵉ + K*η*η̅ + ηᵛᵖ/clock.Δt + ξ^2*H) #ξ*H*∇Q_scalar
        Δλ = Δλ_factor * F_trial

        ##### Test the validity of the return mapping to the smooth portion of the cone
        is_valid = (τ_trial - Gᵛᵉ*Δλ >= 0)

        ##### Use relevant return mapping algorithm

        if is_valid # return to the smooth portion of the cone
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
            A = Δλ_factor
            A_fact = 2Gᵛᵉ * (1.0 - Δλ/(sqrt(2)*norm_ϵᵉdev_t))
            A_factd3 = A_fact/3
            B_fact = 2Gᵛᵉ * (Δλ/(sqrt(2)*norm_ϵᵉdev_t) - Gᵛᵉ*A)
            C_fact = -sqrt(2)*Gᵛᵉ*A*K
            D_fact = K*(1.0 - K*η*η̅*A)
            D = A_fact * Isymdev +
                B_fact * uni_dev ⊗ uni_dev +
                C_fact * (η * uni_dev⊗I2D + η̅ * I2D⊗uni_dev) +
                D_fact * I2D ⊗ I2D #TODO check whether this form of the D term or the one from the book's code is the good one => (D_fact - A_factd3) or just D_fact.

            Δϵᵖ = Δλ*∇Q # plastic strain increment
            Δϵ̅ᵖ = Δλ*ξ # accumulated plastic strain increment
            s.temp_σ = σ_trial - Δλ * Dᵛᵉ⊡∇Q
        else
            α = ξ/η
            β = ξ/η̅
            Δϵᵖ_vol = (p_trial - (C + H*s.ϵ̅ᵖ)*β) / (H*α*β + K)
            Δϵᵖ = 1/3 * Δϵᵖ_vol * I2D
            Δϵ̅ᵖ = α * Δϵᵖ_vol
            s.temp_σ = (p_trial - K*Δϵᵖ_vol)*I2D
            D = K * (1 - K/(K + α*β*H)) * I2D ⊗ I2D
        end
        s.temp_ϵ = ϵ
        s.temp_ϵᵖ = s.ϵᵖ + Δϵᵖ # plastic strain
        s.temp_ϵ̅ᵖ = s.ϵ̅ᵖ + Δϵ̅ᵖ # accumulated plastic strain
        s.temp_σ = σ_trial - Δλ * Dᵛᵉ⊡∇Q
        return s.temp_σ, D
    end
end

loss_of_convexity(u,t,integrator) = (u[1] == 1.0 ? 0 : 1) # must return zero to be triggered

function compute_damaged_stiffness_tensor(r::Rheology,ϵij,D)

    # unpack
    G = r.elasticity.G
    ν = r.elasticity.ν

    # Damage constants
    c1, c2, c3 = compute_c1c2c3(r,D)
    A, B = compute_AB(r,c1,c2,c3)
    A₁ = compute_A1(r,A)
    B₁ = compute_B1(r,B)
    Γ = compute_Γ(r,A₁,B₁)

    # strain invariants
    ϵ = tr(ϵij)
    e = dev(ϵij)
    γ = sqrt(2.0 * e ⊡ e)

    @assert !isnan(G)
    @assert !isnan(c1)
    @assert !isnan(c2)
    @assert !isnan(c3)
    @assert !isnan(A)
    @assert !isnan(B)
    @assert !isnan(A₁)
    @assert !isnan(B₁)
    @assert !isnan(Γ)
    @assert !isnan(ϵ)
    @assert !isnan(γ)
    @assert !isnan((3*(1-2ν))/(2*(1+ν)))
    @assert !isnan(A₁^2/2)
    # println("A₁*B₁*ϵ : ", A₁*B₁*ϵ)
    # println("(2*γ) : ", (2*γ))
    #@assert !isnan(A₁*B₁*ϵ/(2*γ)) # returns NaN, because 0/0

    (γ == 0) && (γ += 1e-9)#zero(typeof(ϵij))) # TODO remove redundancy
    ϵ̂ = ϵij/γ

    # get stiffness factors
    if G/Γ == 0
        Cμ = 0.0
        Cλ = 0.0
        Cσ = 0.0
        Cσσ = 0.0
    else
        Cμ = G/Γ * ( (3*(1-2ν))/(2*(1+ν)) + A₁^2/2 - A₁*B₁*ϵ/(2*γ) )
        Cλ = G/Γ * ( 3*ν/(1+ν) + B₁^2/2 - A₁^2/3 + A₁*B₁*ϵ/γ + 2A₁*B₁*ϵ^3/(9γ^3) )
        Cσ = - G/Γ * ( A₁*B₁ + 2*A₁*B₁*ϵ^2/(3*γ^2) )
        Cσσ = G/Γ * (2A₁*B₁*ϵ/γ)
    end



    if isnan(Cμ)
        println("Cμ = ", Cμ)
        @assert !isnan(G/Γ)
        println("G/Γ = ", G/Γ)
        @assert !isnan( (3*(1-2ν))/(2*(1+ν)) + A₁^2/2 - A₁*B₁*ϵ/(2*γ) )
        println("second term Cμ = ", (3*(1-2ν))/(2*(1+ν)) + A₁^2/2 - A₁*B₁*ϵ/(2*γ) )
        @assert !isnan(G/Γ * ( (3*(1-2ν))/(2*(1+ν)) + A₁^2/2 - A₁*B₁*ϵ/(2*γ) ))
    end
    if isnan(Cλ)
        println("Cλ = ", Cλ)
    end
    @assert !isnan(Cσ)
    @assert !isnan(Cσσ)
    # functional form of the stiffness tensor
    C_func(i,j,k,l) = Cμ * ( δ(k,i)*δ(l,j) + δ(l,i)*δ(k,j) ) +
                      Cλ * ( δ(i,j)*δ(k,l) ) +
                      Cσ * ( ϵ̂[i,j]*δ(k,l) + δ(i,j)*ϵ̂[k,l] ) +
                      Cσσ * (ϵ̂[i,j]*ϵ̂[k,l])

    # assemble the tensor
    return SymmetricTensor{4,3}(C_func)
end

function compute_stress_tangent(ϵij,
                                r::Rheology{T,TD,Nothing,TE,Nothing},
                                s::MaterialState,
                                clock::Clock;
                                noplast = false) where {T,TD<:Damage,TE<:Elasticity}

    # unpack elastic stiffness tensor
    Cᵉ = r.elasticity.Dᵉ

    ##### Evaluation of trial values
    ϵᵉ_trial = ϵij - s.ϵᵖ # trial-elastic_strain
    σ_trial = Cᵉ ⊡ ϵᵉ_trial # trial-stress

    if noplast # during first newton iteration to homogenize strain field
        s.temp_ϵ = ϵij
        s.temp_σ = σ_trial
        return s.temp_σ, Cᵉ
    end

    Δϵ = ϵij - s.ϵ # strain increment
    p_trial = 1/3 * tr(σ_trial) # trial pressure, negative in compression
    s_trial = dev(σ_trial) # trial deviatoric stress
    τ_trial = get_τ(s_trial,r.damage) # effective trial-stress (2nd invariant of the deviatoric stress tensor)
    KI_trial = compute_KI(r,p_trial,τ_trial,s.D) # KI at the end of the timestep if the loading is purely elastic

    if KI_trial < 0.0 # elastic loading over Δt
        s.temp_ϵ = ϵij
        s.temp_σ = σ_trial
        return s.temp_σ, Cᵉ
    else # damage growth over Δt
        # Converged stresses and strains invariants
        ϵ0 = tr(s.ϵ)
        e0 = dev(s.ϵ)
        γ0 = sqrt(2.0 * e0 ⊡ e0)

        σ0 = 1/3 * tr(s.σ)
        s0 = dev(s.σ)
        τ0 = sqrt(0.5 * s0 ⊡ s0)

        # Initial conditions
        u0 = [s.D, ϵ0, γ0, σ0, τ0]
        # Time span
        tspan = (0.0,clock.Δt)
        sub_Δt_ini = tspan[2]/100 # initial subtimestep for RK4 interation algorithm
        damage_reltol = 1e-6
        # strain invariants rates
        dϵdt = tr(Δϵ)/clock.Δt

        e = dev(ϵij)
        γ = sqrt(2.0 * e ⊡ e)
        Δγ = γ - γ0
        dγdt = Δγ/clock.Δt
        # OLD AND WRONG :
        #Δe = dev(Δϵ)
        #dγdt = sqrt(2.0 * Δe ⊡ Δe)/clock.Δt

        # ODE system parameters
        params = [r, dϵdt, dγdt]

        # Hand made RK4
        u_vec, t_vec, cohesion_loss_flag = compute_RK4_adaptative(u0,params,state_system!,tspan,sub_Δt_ini,damage_reltol)

        # TEST CODE
        #println("t_vec : ", t_vec)
        #println("D_vec : ", [u[1] for u in u_vec])
        #lineplot(t_vec,u_vec[1])


        # DiffEq solve alternatively
        # # Problem instantiation
        # state_problem = ODEProblem(state_system!,u0,tspan,params)
        #
        # # the ODE system is only solved until D reaches 1.0 => loss of convexity of the free energy function.
        # termination_callback = ContinuousCallback(loss_of_convexity,terminate!)
        # sol = OrdinaryDiffEq.solve(state_problem, Tsit5(), callback=termination_callback)

        # TODO for damaged-plastic material, the idea is to check if the integration has been done over the whole Δt. If it stopped before, it means loss of cohesion of the material. A plastic update is therefore required in that case.


        # update states
        s.temp_D = u_vec[end][1]#sol[1,end]

        if s.temp_D > s.D
            s.temp_σ = compute_σij(r,s.temp_D,ϵij) # changed from ϵᵉ_trial to ϵ

            # TEST
            σ_temp = 1/3 * tr(s.temp_σ) # trial pressure, negative in compression
            s_temp = dev(s.temp_σ) # trial deviatoric stress
            τ_temp = get_τ(s_temp,r.damage)

            #println("sigma diff : ", σ_temp - u_vec[end][4])
            #println("tau diff : ", τ_temp - u_vec[end][5])
            invCᵉ = get_elastic_compliance_tensor(r)
            Δϵᵖ = Δϵ - invCᵉ⊡(s.temp_σ - s.σ) # damage strain increment
            s.temp_ϵ = ϵij
            s.temp_ϵᵖ = s.ϵᵖ + Δϵᵖ # plastic strain
            # Construction of the tangent operator ∂σ_n/∂ϵ_n at {D = D_n+1}:
            C = compute_damaged_stiffness_tensor(r,s.ϵ,s.D) #TODO s.D or S.temp_D
            return s.temp_σ, C
        else
            s.temp_ϵ = ϵij
            s.temp_σ = σ_trial
            return s.temp_σ, Cᵉ
        end
        # EQUALITY TEST TODO investigate that
        #Δσ = C ⊡ Δϵ
        # println("Δϵ : ", Δϵ)
        # println("tangent : ", C)
        # println("Δσ from tangent : ", Δσ)
        # println("Δσ from constitutive relationships : ", s.temp_σ - s.σ)
        #Δσ = s.temp_σ - s.σ
        #println("Δσ : ", Δσ)
        #println("Δϵ : ", Δϵ)
        #C2 = Δσ ⊗ pinv(Δϵ)
        #Δσ_C2 = C2 ⊡ Δϵ
        #println("Δσ const rel : ", Δσ)
        #println("Δσ from Δσ ⊗ inv(Δϵ) : ", Δσ_C2)


        #Δϵ̅ᵖ = Δλ*ξ#Δλ*∇Q_scalar # accumulated plastic strain increment

        #s.temp_ϵ̅ᵖ = s.ϵ̅ᵖ + Δϵ̅ᵖ # accumulated plastic strain
        #s.temp_σ = σ_trial - Δλ * Dᵉ⊡∇Q
    end
end

function compute_stress_tangent(ϵij,
                                r::Rheology{T,TD,Nothing,TE,TP},
                                s::MaterialState,
                                clock::Clock;
                                noplast = false) where {T,TD<:Damage,TE<:Elasticity,TP}

    # unpack elastic stiffness tensor
    Cᵉ = r.elasticity.Dᵉ

    ##### Evaluation of trial values
    ϵᵉ_trial = ϵij - s.ϵᵖ # trial-elastic_strain
    σ_trial = Cᵉ ⊡ ϵᵉ_trial # trial-stress

    if noplast # during first newton iteration to homogenize strain field
        s.temp_ϵ = ϵij
        s.temp_σ = σ_trial
        return s.temp_σ, Cᵉ
    end

    Δϵ = ϵij - s.ϵ # strain increment
    p_trial = 1/3 * tr(σ_trial) # trial pressure, negative in compression
    s_trial = dev(σ_trial) # trial deviatoric stress
    τ_trial = get_τ(s_trial,r.damage) # effetive trial-stress (2nd invariant of the deviatoric stress tensor)
    KI_trial = compute_KI(r,p_trial,τ_trial,s.D) # KI at the end of the timestep if the loading is purely elastic

    if KI_trial < 0.0 # elastic loading over Δt
        s.temp_ϵ = ϵij
        s.temp_σ = σ_trial
        return s.temp_σ, Cᵉ
    else # damage growth over Δt
        # Converged stresses and strains invariants
        ϵ0 = tr(s.ϵ)
        e0 = dev(s.ϵ)
        γ0 = sqrt(2.0 * e0 ⊡ e0)

        σ0 = 1/3 * tr(s.σ)
        s0 = dev(s.σ)
        τ0 = sqrt(0.5 * s0 ⊡ s0)

        # Initial conditions
        u0 = [s.D, ϵ0, γ0, σ0, τ0]
        # Time span
        tspan = (0.0,clock.Δt)
        sub_Δt_ini = tspan[2]/100 # initial subtimestep for RK4 interation algorithm
        damage_reltol = 1e-6
        # strain invariants rates
        dϵdt = tr(Δϵ)/clock.Δt

        e = dev(ϵij)
        γ = sqrt(2.0 * e ⊡ e)
        Δγ = γ - γ0
        dγdt = Δγ/clock.Δt
        #Δe = dev(Δϵ)
        #dγdt = sqrt(2.0 * Δe ⊡ Δe)/clock.Δt

        # ODE system parameters
        params = [r, dϵdt, dγdt]

        # Hand made RK4
        u_vec, t_vec, cohesion_loss_flag = compute_RK4_adaptative(u0,params,state_system!,tspan,sub_Δt_ini,damage_reltol)

        # TEST CODE
        #println("t_vec : ", t_vec)
        #println("D_vec : ", [u[1] for u in u_vec])
        #lineplot(t_vec,u_vec[1])


        # DiffEq solve alternatively
        # # Problem instantiation
        # state_problem = ODEProblem(state_system!,u0,tspan,params)
        #
        # # the ODE system is only solved until D reaches 1.0 => loss of convexity of the free energy function.
        # termination_callback = ContinuousCallback(loss_of_convexity,terminate!)
        # sol = OrdinaryDiffEq.solve(state_problem, Tsit5(), callback=termination_callback)

        # TODO for damaged-plastic material, the idea is to check if the integration has been done over the whole Δt. If it stopped before, it means loss of cohesion of the material. A plastic update is therefore required in that case.

        # Handle loss of cohesion => branching to a plastic behavior
        if (TP <: Plasticity) & cohesion_loss_flag
            r_nodamage = Rheology(damage = Nothing,
                                  viscosity = r.viscosity,
                                  elasticity = r.elasticity,
                                  plasticity = r.plasticity)
                                  # TODO tester avec ou sans reset de plastic strain
            compute_stress_tangent(ϵij, r_nodamage, s, clock)
        end

        # update state
        s.temp_D = u_vec[end][1]#sol[1,end]

        if s.temp_D > s.D
            s.temp_σ = compute_σij(r,s.temp_D,ϵij) # changed from ϵᵉ_trial to ϵ

            ## TEST ##
            if s.temp_D - s.D > 0.0001
                p_temp = 1/3 * tr(s.temp_σ) # pressure, negative in compression
                s_temp = dev(s.temp_σ) # deviatoric stress
                τ_temp = get_τ(s_temp,r.damage)
                println("p diff : ", p_temp - u_vec[end][4])
                println("tau diff : ", τ_temp - u_vec[end][5])
            end
            invCᵉ = get_elastic_compliance_tensor(r)
            Δϵᵖ = Δϵ - invCᵉ⊡(s.temp_σ - s.σ) # damage strain increment
            s.temp_ϵ = ϵij
            s.temp_ϵᵖ = s.ϵᵖ + Δϵᵖ # plastic strain
            # Construction of the tangent operator ∂σ_n/∂ϵ_n at {D = D_n+1}:
            C = compute_damaged_stiffness_tensor(r,s.ϵ,s.D) #TODO check s.D or temp_D
            return s.temp_σ, C
        else
            s.temp_ϵ = ϵij
            s.temp_σ = σ_trial
            return s.temp_σ, Cᵉ
        end


        # EQUALITY TEST TODO investigate that
        #Δσ = C ⊡ Δϵ
        # println("Δϵ : ", Δϵ)
        # println("tangent : ", C)
        # println("Δσ from tangent : ", Δσ)
        # println("Δσ from constitutive relationships : ", s.temp_σ - s.σ)
        #Δσ = s.temp_σ - s.σ
        #println("Δσ : ", Δσ)
        #println("Δϵ : ", Δϵ)
        #C2 = Δσ ⊗ pinv(Δϵ)
        #Δσ_C2 = C2 ⊡ Δϵ
        #println("Δσ const rel : ", Δσ)
        #println("Δσ from Δσ ⊗ inv(Δϵ) : ", Δσ_C2)


        #Δϵ̅ᵖ = Δλ*ξ#Δλ*∇Q_scalar # accumulated plastic strain increment

        #s.temp_ϵ̅ᵖ = s.ϵ̅ᵖ + Δϵ̅ᵖ # accumulated plastic strain
        #s.temp_σ = σ_trial - Δλ * Dᵉ⊡∇Q
    end
end

function compute_stress_tangent(ϵij,
                                r::Rheology{T,TD,TV,TE,TP},
                                s::MaterialState,
                                clock::Clock;
                                noplast = false) where {T,TD<:Damage,TV<:Viscosity,TE<:Elasticity,TP<:Plasticity}

    ##### viscoelastic tangent operator
    Gᵛᵉ = 1/(1/r.elasticity.G + clock.Δt/r.viscosity.η)
    λᵛᵉ = λ_from_KG(r.elasticity.K,Gᵛᵉ)
    Cᵛᵉ = get_elastic_stiffness_tensor(Gᵛᵉ,λᵛᵉ)

    # trial stress
    Gratio = Gᵛᵉ/r.elasticity.G
    I2D = SymmetricTensor{2,3}(δ) # second order identity tensor
    σ⁰ = (1/3)*tr(s.σ)*I2D + Gratio*dev(s.σ)
    ϵᵉ_trial = ϵij - s.ϵᵖ # trial-elastic_strain
    Δϵ = ϵij - s.ϵ # strain increment
    σ_trial = σ⁰ + Cᵛᵉ ⊡ (ϵij - s.ϵ) # trial stress

    ##### Evaluation of trial values OLD
    # ϵᵉ_trial = ϵ - s.ϵᵖ # trial-elastic_strain
    # σ_trial = Cᵉ ⊡ ϵᵉ_trial # trial-stress

    if noplast # during first newton iteration to homogenize strain field
        s.temp_ϵ = ϵij
        s.temp_σ = σ_trial
        return s.temp_σ, Cᵛᵉ
    end

    p_trial = 1/3 * tr(σ_trial) # trial pressure, negative in compression
    s_trial = dev(σ_trial) # trial deviatoric stress
    τ_trial = get_τ(s_trial,r.damage) # effetive trial-stress (2nd invariant of the deviatoric stress tensor)
    KI_trial = compute_KI(r,p_trial,τ_trial,s.D) # KI at the end of the timestep if the loading is purely elastic

    if KI_trial < 0.0 # elastic loading over Δt
        s.temp_ϵ = ϵij
        s.temp_σ = σ_trial
        return s.temp_σ, Cᵛᵉ
    else # damage growth over Δt
        # Converged stresses and strains invariants
        ϵ0 = tr(s.ϵ)
        e0 = dev(s.ϵ)
        γ0 = sqrt(2.0 * e0 ⊡ e0)

        σ0 = 1/3 * tr(s.σ)
        s0 = dev(s.σ)
        τ0 = sqrt(0.5 * s0 ⊡ s0)

        # Initial conditions
        u0 = [s.D, ϵ0, γ0, σ0, τ0]
        # Time span
        tspan = (0.0,clock.Δt)
        sub_Δt_ini = tspan[2]/100 # initial subtimestep for RK4 interation algorithm
        damage_reltol = 1e-6
        # strain invariants rates
        dϵdt = tr(Δϵ)/clock.Δt
        e = dev(ϵij)
        γ = sqrt(2.0 * e ⊡ e)
        Δγ = γ - γ0
        dγdt = Δγ/clock.Δt
        # OLD AND WRONG :
        #Δe = dev(Δϵ)
        #dγdt = sqrt(2.0 * Δe ⊡ Δe)/clock.Δt

        # ODE system parameters
        params = [r, dϵdt, dγdt]

        # Hand made RK4
        u_vec, t_vec, cohesion_loss_flag = compute_RK4_adaptative(u0,params,state_system!,tspan,sub_Δt_ini,damage_reltol)

        # TEST CODE
        #println("t_vec : ", t_vec)
        #println("D_vec : ", [u[1] for u in u_vec])
        #lineplot(t_vec,u_vec[1])


        # DiffEq solve alternatively
        # # Problem instantiation
        # state_problem = ODEProblem(state_system!,u0,tspan,params)
        #
        # # the ODE system is only solved until D reaches 1.0 => loss of convexity of the free energy function.
        # termination_callback = ContinuousCallback(loss_of_convexity,terminate!)
        # sol = OrdinaryDiffEq.solve(state_problem, Tsit5(), callback=termination_callback)

        # Handle loss of cohesion => branching to a plastic behavior
        if u_vec[end][1] >= 1
            r_nodamage = Rheology(damage = Nothing,
                                  viscosity = r.viscosity,
                                  elasticity = r.elasticity,
                                  plasticity = r.plasticity)
            compute_stress_tangent(ϵij, r_nodamage,s,clock ; noplast = false)
        end

        # update temporary state
        s.temp_D = u_vec[end][1]

        if s.temp_D > s.D
            s.temp_σ = compute_σij(r,s.temp_D,ϵij) # changed from ϵᵉ_trial to ϵ
            #invCᵉ = get_elastic_compliance_tensor(r) # TODO add a viscoelastic compliance tensor getter
            Δϵᵖ = Δϵ - inv(Cᵛᵉ)⊡(s.temp_σ - s.σ) # damage strain increment
            s.temp_ϵ = ϵij
            s.temp_ϵᵖ = s.ϵᵖ + Δϵᵖ # plastic strain
            # Construction of the tangent operator ∂σ_n/∂ϵ_n at {D = D_n+1}:
            C = compute_damaged_stiffness_tensor(r,s.ϵ,s.temp_D)
            return s.temp_σ, C
        else
            s.temp_ϵ = ϵij
            s.temp_σ = σ_trial
            return s.temp_σ, Cᵛᵉ
        end
        # EQUALITY TEST TODO investigate that
        #Δσ = C ⊡ Δϵ
        # println("Δϵ : ", Δϵ)
        # println("tangent : ", C)
        # println("Δσ from tangent : ", Δσ)
        # println("Δσ from constitutive relationships : ", s.temp_σ - s.σ)
        #Δσ = s.temp_σ - s.σ
        #println("Δσ : ", Δσ)
        #println("Δϵ : ", Δϵ)
        #C2 = Δσ ⊗ pinv(Δϵ)
        #Δσ_C2 = C2 ⊡ Δϵ
        #println("Δσ const rel : ", Δσ)
        #println("Δσ from Δσ ⊗ inv(Δϵ) : ", Δσ_C2)


        #Δϵ̅ᵖ = Δλ*ξ#Δλ*∇Q_scalar # accumulated plastic strain increment

        #s.temp_ϵ̅ᵖ = s.ϵ̅ᵖ + Δϵ̅ᵖ # accumulated plastic strain
        #s.temp_σ = σ_trial - Δλ * Dᵉ⊡∇Q
    end
end

function compute_jacobian(f!, u ; eps=1e-9)
    Iv = Int[]
    Jv = Int[]
    Vv = Float64[]
    len_u = length(u)
    δu = Vector{Float64}(undef,len_u)
    half_δu = similar(δu)
    ∂f∂u_j = similar(δu)
    f_up = similar(δu)
    f_down = similar(δu)
    @inbounds for j in 1:len_u
        fill!(δu,0.0)
        δu[j] = eps
        half_δu .= δu./2
        f!(f_up, u .+ half_δu)
        f!(f_down, u .- half_δu)
        @. ∂f∂u_j = ( f_up - f_down ) /  δu[j]
        for i in 1:len_u
            if ∂f∂u_j[i] != 0.0
                push!(Iv,i)
                push!(Jv,j)
                push!(Vv,∂f∂u_j[i])
            end
        end
    end
    J = sparse(Iv,Jv,Vv,len_u,len_u)
    return J
end
