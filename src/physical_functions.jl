### PHYSICAL FUNCTIONS ###
include("RK4_functions.jl")

get_elastic_stiffness_tensor(G,λ) = SymmetricTensor{4, 3}( (i,j,k,l) -> Dᵉ_func(i,j,k,l,G,λ))
get_elastic_stiffness_tensor(e::Elasticity) = SymmetricTensor{4, 3}( (i,j,k,l) -> Dᵉ_func(i,j,k,l,e.G,e.λ))
get_elastic_compliance_tensor(e::Elasticity) = SymmetricTensor{4, 3}( (i,j,k,l) -> Cᵉ_func(i,j,k,l,e.G,e.λ))

get_elastic_stiffness_tensor(r::Rheology) = get_elastic_stiffness_tensor(r.elasticity)
get_elastic_compliance_tensor(r::Rheology) = get_elastic_compliance_tensor(r.elasticity)
# Second invariant of the deviatoric stress tensor :
get_τ(s ; plas = :DruckerPrager) = sqrt(0.5 * s ⊡ s)
get_τ(s,::DruckerPrager) = sqrt(0.5 * s ⊡ s)
get_τ(s,::Damage) = sqrt(0.5 * s ⊡ s)
get_τ(s,::VonMises) = sqrt(3/2 * s ⊡ s)

get_τ(s,r::Rheology) = get_τ(s::AbstractTensor,r.plasticity)

### Damage functions ###

free_energy_convexity(r,D,A₁,B₁) = 1/compute_Γ(r,A₁,B₁) > 0 ? true : false

# eq 16 Bhat2012 & 2016 & notes (because in Bhat2011 c2 isn't the same form as in Harsha's notes) :
#compute_c1(d::Damage,D) = sqrt(1-cos(d.ψ)^2)/(π*cos(d.ψ)^(3/2)*((D/d.D₀)^(1/3) - 1 + d.β/cos(d.ψ))^(3/2))
function compute_c1(d::Damage,D)
    α = cosd(d.ψ)
    @assert α > 0
    if d.D₀>0
        @assert (D/d.D₀) >= 1
        sqrt(1-α^2)/(π*α^(3/2)*((D/d.D₀)^(1/3) - 1 + d.β/α)^(3/2))
    elseif (d.D₀==0) & (D==0)
        return sqrt(1-α^2)/(π*α^(3/2)*(d.β/α)^(3/2)) # lim D->0 D/D0 = 1
      end
end
# Perol&Bhat2016 : 1/α  or  Harsha's notes : 1/α^2 ???
function compute_c2(d::Damage,D)
    α = cosd(d.ψ)
    (d.D₀==0) && (return 0.0)
    return (sqrt(1 - α^2)/α^2) * (d.D₀^(2/3)/(1 - D^(2/3))) ###
  end

function compute_c3(d::Damage,D)
    α = cosd(d.ψ)
    @assert α > 0
    (d.D₀ == D) && (return 0.0)
    @assert (D/d.D₀) >= 1
    (2sqrt(α)/π)*((D/d.D₀)^(1/3) - 1)^(1/2)
end

compute_c1(r::Rheology,D) = compute_c1(r.damage::Damage,D)
compute_c2(r::Rheology,D) = compute_c2(r.damage::Damage,D)
compute_c3(r::Rheology,D) = compute_c3(r.damage::Damage,D)

function compute_c1c2c3(d::Damage,D)
    c1 = compute_c1(d,D)
    c2 = compute_c2(d,D)
    c3 = compute_c3(d,D)
    return c1, c2, c3
end

compute_c1c2c3(r::Rheology,D) = compute_c1c2c3(r.damage::Damage,D)

# eq 15 Bhat2012 (A1 : *c2*c3), Perol&Bhat2016 (A1 : ...*c2)*c3):
# Perol&Bhat2016 is the corrected version, and the one implemented
compute_A(r::Rheology,c1,c2,c3) = r.damage.μ*c1 + (1.0 + r.damage.μ*c2)*c3
compute_A(d::Damage,c1,c2,c3) = d.μ*c1 + (1.0 + d.μ*c2)*c3
compute_B(c1,c2,c3) = c1 + c2*c3

function compute_AB(d::Damage,c1,c2,c3)
    A = compute_A(d,c1,c2,c3)
    B = compute_B(c1,c2,c3)
    return A, B
end
compute_AB(r::Rheology,c1,c2,c3) = compute_AB(r.damage,c1,c2,c3)
compute_AB(r::Rheology,D) = compute_AB(r,compute_c1c2c3(r,D)...)


# eq 11 in Harsha's notes :
compute_A1(r::Rheology,A) = A * sqrt((π*r.damage.D₀*(1 - r.elasticity.ν))/cosd(r.damage.ψ)^3)
compute_B1(r::Rheology,B) = B * sqrt((π*r.damage.D₀*(1 - r.elasticity.ν))/cosd(r.damage.ψ)^3)

function compute_A1B1(r::Rheology,A,B)
    A1 = compute_A1(r,A)
    B1 = compute_B1(r,B)
    return A1, B1
end

compute_A1B1(r::Rheology,D) = compute_A1B1(r,compute_AB(r,D)...)

function compute_a1(B1,Γ)
    return (1/Γ)*(1 + B1^2/2)
end

function compute_b1(A1,B1,Γ)
    return -(1/Γ)*((A1*B1)/2)
end

function compute_b2(r::Rheology,A1,Γ)
    return (1/Γ)*(A1^2/2 + (3*(1-2r.elasticity.ν))/(2*(1+r.elasticity.ν)))
end

function compute_Γ(r::Rheology,A₁,B₁)
    ν = r.elasticity.ν
    return (3*(1-2ν))/(2*(1+ν)) + (3*(1-2ν)*B₁^2)/(4*(1+ν)) + A₁^2/2
end

function compute_Γ(r::Rheology,D)
    c1 = compute_c1(r,D)
    c2 = compute_c2(r,D)
    c3 = compute_c3(r,D)
    A = compute_A(r,c1,c2,c3)
    B = compute_B(c1,c2,c3)
    A₁ = compute_A1(r,A)
    B₁ = compute_B1(r,B)
    ν = r.elasticity.ν
    return (3*(1-2ν))/(2*(1+ν)) + (3*(1-2ν)*B₁^2)/(4*(1+ν)) + A₁^2/2
end

function compute_dc1dD(r::Rheology,D)
    d = r.damage
    α = cosd(d.ψ)
    term1 = (-sqrt(1-α^2))/(2*π*α^(3/2)*D^(2/3)*d.D₀^(1/3))
    term2 = ((D/d.D₀)^(1/3) - 1 + (d.β/α))^(-5/2)
    return term1 * term2
end

function compute_dc2dD(r::Rheology,D)
    d = r.damage
    α = cosd(d.ψ)
    term1 = (2*sqrt(1-α^2)*d.D₀^(2/3))/(3*α^2*D^(1/3))
    term2 = (1 - D^(2/3))^(-2)
    return term1 * term2
end

function compute_dc3dD(r::Rheology,D)
    d = r.damage
    α = cosd(d.ψ)
    term1 = (sqrt(α))/(3*π*D^(2/3)*d.D₀^(1/3))
    term2 = ((D/d.D₀)^(1/3) - 1)^(-1/2)

    # ensure term2 isn't Infinity
    (term2 == Inf) && (term2 = 1e15)
    return term1 * term2
end

function compute_dA1dD(r::Rheology,dc1dD,dc2dD,dc3dD,c2,c3)
    d = r.damage
    α = cosd(d.ψ)
    term1 = sqrt((π*d.D₀*(1 - r.elasticity.ν))/(α^3))
    term2 = d.μ*dc1dD + dc3dD + d.μ*c2*dc3dD + d.μ*c3*dc2dD
    return term1 * term2
end

function compute_dB1dD(r::Rheology,dc1dD,dc2dD,dc3dD,c2,c3)
    d = r.damage
    α = cosd(d.ψ)
    term1 = sqrt((π*d.D₀*(1 - r.elasticity.ν))/(α^3))
    term2 = dc1dD + c2*dc3dD + c3*dc2dD
    return term1 * term2
end

function compute_dΓdD(r::Rheology,A1,B1,dA1dD,dB1dD)
    return ((3*(1-2r.elasticity.ν))/(2*(1+r.elasticity.ν)))*B1*dB1dD + A1*dA1dD
end

function compute_dA1dt(r::Rheology,dc1dt,dc2dt,dc3dt,c2,c3)
    d = r.damage
    α = cosd(d.ψ)
    term1 = sqrt((π*d.D₀*(1 - r.elasticity.ν))/(α^3))
    term2 = d.μ*dc1dt + dc3dt + d.μ*c2*dc3dt + d.μ*c3*dc2dt
    return term1 * term2
end

function compute_dB1dt(r::Rheology,dc1dt,dc2dt,dc3dt,c2,c3)
    d = r.damage
    α = cosd(d.ψ)
    term1 = sqrt((π*d.D₀*(1 - r.elasticity.ν))/(α^3))
    term2 = dc1dt + c2*dc3dt + c3*dc2dt
    return term1 * term2
end

function compute_dΓdt(r::Rheology,A1,B1,dA1dt,dB1dt)
    return ((3*(1-2r.elasticity.ν))/(2*(1+r.elasticity.ν)))*B1*dB1dt + A1*dA1dt
end

# function compute_dΓdt(p::RockParams,A1,B1,dA1dt,dB1dt)
#     return ((3*(1-2p.ν))/(2*(1+p.ν)))*B1*dB1dt + A1*dA1dt
# end

function compute_da1dt(B1,Γ,dB1dt,dΓdt)
    return -(dΓdt/Γ^2)*(1 + B1^2/2) + (B1*dB1dt)/Γ
end

function compute_db1dt(A1,B1,Γ,dA1dt,dB1dt,dΓdt)
    return (dΓdt/Γ^2)*((A1*B1)/2) - (1/2Γ)*(dA1dt*B1 + A1*dB1dt)
end

function compute_db2dt(r::Rheology,A1,Γ,dA1dt,dΓdt)
    return -(dΓdt/Γ^2)*(A1^2/2 + (3*(1-2r.elasticity.ν))/(2*(1+r.elasticity.ν))) + (2/Γ)*A1*dA1dt
end

function compute_dσdt(r::Rheology,a1,b1,da1dt,db1dt,ϵ,γ,dϵdt,dγdt)
    return r.elasticity.G * (da1dt*ϵ + a1*dϵdt + db1dt*γ + b1*dγdt)
end

function compute_dτdt(r::Rheology,b1,b2,db1dt,db2dt,ϵ,γ,dϵdt,dγdt)
    return r.elasticity.G * (db1dt*ϵ + b1*dϵdt + db2dt*γ + b2*dγdt)
end

# Damage growth and KI

function compute_KI(r::Rheology,σ,τ,A,B)
    return (A*σ + B*τ) * sqrt(π*r.damage.a)
end

function compute_KI(r::Rheology,σ,τ,D)
    c1, c2, c3 = compute_c1c2c3(r,D)
    A, B = compute_AB(r,c1,c2,c3)
    # println("c1 : ",c1)
    # println("c2 : ",c2)
    # println("c3 : ",c3)
    # println("σ : ",σ)
    # println("τ : ",τ)
    return (A*σ + B*τ) * sqrt(π*r.damage.a)
end

function compute_KI(d::Damage,σij::SymmetricTensor,D)
    c1, c2, c3 = compute_c1c2c3(d,D)
    A, B = compute_AB(d,c1,c2,c3)
    p = 1/3 * tr(σij) # trial pressure, negative in compression
    sij = dev(σij) # trial deviatoric stress
    τ = get_τ(sij,d)
    return (A*p + B*τ) * sqrt(π*d.a)
end
compute_KI(r::Rheology,σij,D) = compute_KI(r.damage,σij,D)

function compute_dDdl(r::Rheology,D)
    d = r.damage
    return (3*D^(2/3)*d.D₀^(1/3))/(cosd(d.ψ)*d.a)
end

function compute_subcrit_damage_rate(r::Rheology, KI, D)
    ((KI <= 0) | (D >= 1)) && (return 0.0)
    d = r.damage
    e = r.elasticity
    ρ = 2700 ##### TODO better
    Vs = sqrt(e.G/ρ)
    Vr = Vs * (0.862 + 1.14e.ν)/(1 + e.ν)

    dDdl = compute_dDdl(r,D) # damage derivative wrt crack length
    dldt = min(d.l̇₀*(KI/d.K₁c)^(d.n),Vr)  # cracks growth rate
    @assert dDdl * dldt >= 0
    return dDdl * dldt
end

function compute_subcrit_damage_rate(r::Rheology, σij::SymmetricTensor, D)
    KI = compute_KI(r,σij,D)
    ((KI <= 0) | (D >= 1)) && (return 0.0)
    d = r.damage
    e = r.elasticity
    ρ = 2700 ##### TODO better
    Vs = sqrt(e.G/ρ)
    Vr = Vs * (0.862 + 1.14e.ν)/(1 + e.ν)

    dDdl = compute_dDdl(r,D) # damage derivative wrt crack length
    dldt = min(d.l̇₀*(KI/d.K₁c)^(d.n),Vr)  # cracks growth rate
    @assert dDdl * dldt >= 0
    return dDdl * dldt
end

# dDdt_thresh=1e-10
# KIC_max=1e8
# ΔKI = d.K₁c*(dDdt_thresh/(1.5e3*0.24))^(1/d.n)
# KI = -1e6:100:3e6
# x = 2.5.*KI./ΔKI
# KIC = (1.0 .- tanh.(x)).*(KIC_max - d.K₁c).*0.5 .+ d.K₁c

function compute_subcrit_damage_rate_variable_KIC(r::Rheology, KI, D ; dDdt_thresh=1e-40, KIC_max=1e15 )
    #((KI <= 0) | (D >= 1)) && (return 0.0)
    d = r.damage
    e = r.elasticity
    ρ = 2700 ##### TODO better
    Vs = sqrt(e.G/ρ)
    Vr = Vs * (0.862 + 1.14e.ν)/(1 + e.ν)

    dDdl = compute_dDdl(r,D) # damage derivative wrt crack length

    # get KI length scale of transition from high KIC to real KIC
    #println("KI : ",KI)
    KIC = compute_KIC(r, KI, dDdt_thresh, KIC_max)

    dldt = min(d.l̇₀*(abs(KI)/KIC)^(d.n),Vr)  # cracks growth rate

    try @assert dDdl * dldt >= 0
    catch e
        println("KI : ", KI)
        println("KIC : ", KIC)
        throw(e)
    end

    return dDdl * dldt
end

function compute_KIC(r, KI, dDdt_thresh, KIC_max)
    d = r.damage
    ΔKI = d.K₁c*(dDdt_thresh/(1.5e3*d.l̇₀))^(1/d.n)
    x = 2.5*KI/ΔKI
    #println("x : ", x)
    #println("tanh(x) : ", tanh(x))
    KIC = (1.0 - tanh(x))*(KIC_max - d.K₁c)*0.5 + d.K₁c
    return KIC
end

function compute_subcrit_damage_rate(r::Rheology, σ, τ, D)
    KI = compute_KI(r,σ,τ,D)
    # return zero if KI is negative
    ((KI <= 0) | (D >= 1)) && (return 0.0)

    d = r.damage
    e = r.elasticity
    ρ = 2700 ##### TODO better
    Vs = sqrt(e.G/ρ)
    Vr = Vs * (0.862 + 1.14e.ν)/(1 + e.ν)

    dDdl = compute_dDdl(r,D) # damage derivative wrt crack length
    dldt = min(d.l̇₀*(KI/d.K₁c)^(d.n),Vr)  # cracks growth rate  speed
    # println("d.l̇₀ = ",d.l̇₀)
    # println("d.K₁c = ",d.K₁c)
    # println("d.n = ",d.n)
    # println("Vr = ",Vr)
    # println("dldt = ",dldt)
    @assert dDdl >= 0

    return dDdl * dldt
end

function get_damage_constrained_Δt(model,u,ΔD_max)
    nu,nD = getnbasefunctions.(model.cellvalues_tuple)
    nqp = getnquadpoints(model.cellvalues_tuple[2])
    Δt_max = 1e9
    dDdt_max = 0.0
    for cell in CellIterator(model.dofhandler)
        cellid = cell.current_cellid.x
        states = model.material_state[cellid]
        r = model.material_properties[cellid]

        cell_dofs = celldofs(cell)
        ue = u[cell_dofs]
        De = exp.(ue[nu+1:end])
        for qp in nqp
            state = states[qp]
            D = function_value(model.cellvalues_tuple[2],qp,De)
            KI = compute_KI(r,state.temp_σ,D)
            dDdt = compute_subcrit_damage_rate(r, KI, D)
            Δt_max = min(ΔD_max/dDdt,Δt_max)
            dDdt_max = max(dDdt,dDdt_max)
        end
    end
    #println("dDdt_max based on converged σ and D = ",dDdt_max)
    return Δt_max
end

# function get_damage_constrained_Δt(model,ΔD_max)
#     nqp = getnquadpoints(model.cellvalues_tuple[1])
#     KI = 0.0
#     for cell in CellIterator(model.dofhandler)
#         cellid = cell.current_cellid.x
#         r = model.material_properties[cellid]
#         for qp in nqp
#             state = model.material_state[cellid][qp]
#             KI = max(KI,compute_KI(r,state.temp_σ,D))
#         end
#     end
#     KI_max =
#     D = function_value(model.cellvalues_tuple[2],qp,De)
#     KI = compute_KI(r,state.temp_σ,D)
#     dDdt = compute_subcrit_damage_rate(r, KI, D)
#     Δt_max = min(ΔD_max/dDdt,Δt_max)
#     dDdt_max = max(dDdt,dDdt_max)
#     println("dDdt_max based on converged σ and D = ",dDdt_max)
#     return Δt_max
# end
# function compute_subcrit_damage_rate(r::Rheology, σ::T, D) where {T<:AbstractArray}
#     d = r.damage
#     e = r.elasticity
#     G = e.E / (2*(1 + e.ν))
#     Vs = sqrt(G/2700) ##### TODO better
#
#     KI = compute_KI(r,σ,D)
#     dDdl = compute_dDdl(r,D) # damage derivative wrt crack length
#     dldt = min(d.l̇₀*(KI/d.K₁c)^(d.n),Vs)  # cracks growth rate #TODO should be rayleigh wave speed
#     return dDdl * dldt
# end


# function compute_dσijdt(r,ϵij, D, dϵijdt, dDdt)
#
#     # strain decomposition
#     ϵ = tr(ϵij)
#     e = dev(ϵij)
#     γ = sqrt(2.0 * e ⊡ e)
#     dϵdt = tr(dϵijdt)
#     dedt = dev(dϵijdt)
#     dγdt = sqrt(2.0 * dedt ⊡ dedt)
#
#     # damage constants
#     c1, c2, c3 = compute_c1c2c3(r,D)
#     A, B = compute_AB(r,c1,c2,c3)
#     A1 = compute_A1(r,A)
#     B1 = compute_B1(r,B)
#     Γ = compute_Γ(r,A1,B1)
#
#     # derivatives
#     dc1dD = compute_dc1dD(r,D)
#     dc2dD = compute_dc2dD(r,D)
#     dc3dD = compute_dc3dD(r,D)
#     dA1dD = compute_dA1dD(r,dc1dD,dc2dD,dc3dD,c2,c3)
#     dB1dD = compute_dB1dD(r,dc1dD,dc2dD,dc3dD,c2,c3)
#     dΓdD = compute_dΓdD(r,A₁,B₁,dA1dD,dB1dD)
#
#     dA1dt = dA1dD * dDdt
#     dB1dt = dB1dD * dDdt
#     dΓdt = dΓdD * dDdt
#     # unpack
#     G = r.elasticity.G
#     ν = r.elasticity.ν
#
#     # Identity matrix
#     Id = SymmetricTensor{2,3}(δ)
#     # dσijdt, with σij = F*H
#     F = G / Γ
#     H = ( (3*(1-2ν))/(1+ν) + A1^2 - A1*B1*ϵ/γ )*ϵij +
#         ( (3ν/(1+ν) + B1^2/2 - A1^2/3 + A1*B1*ϵ/(3γ))*ϵ - A1*B1*γ/2 ) * Id
#
#     dFdt = -G/Γ^2 * dΓdt
#     dHdt = ( 2A1*dA1dt - (dA1dt*B1 + A1*dB1dt)*ϵ/γ - A1*B1*(dϵdt/γ - ϵ*dγdt/γ^2) )*ϵij +
#            ( (3*(1-2ν))/(1+ν) + A1^2 - A1*B1*ϵ/γ )*dϵijdt +
#            ( (B1*dB1dt - (2/3)*A1*dA1dt + (1/3)*((dA1dt*B1 + A1*dB1dt)*ϵ/γ - A1*B1*(dϵdt/γ - ϵ*dγdt/γ^2)))*ϵ +
#              (3ν/(1+ν) + B1^2/2 - A1^2/3 + A1*B1*ϵ/(3γ))*dϵdt -
#              0.5*((dA1dt*B1 + A1*dB1dt)*γ + A1*B1*dγdt) ) * Id
#
#     return dFdt*H + F*dHdt
# end

# function state_tensors_system!(du,u,p,t)
#     D, ϵ, σ = u
#     rheology, Δt, Δϵ = p
#     KI = compute_KI(rheology,σ,D)
#     du[1] = dD = (KI > 0) ? compute_subcrit_damage_rate(rheology, σ, D) : 0.0
#     du[2] = dϵ = Δϵ/Δt
#     du[3] = dσ = compute_dσijdt(rheology, ϵ, D, du[2], du[1])
# end
function compute_ϵ_oop(r::Rheology,A1,B1,σᵢⱼ::SymmetricTensor)

    σoop = σᵢⱼ[3,3]
    σ = 1/3*tr(σᵢⱼ)
    sᵢⱼ = dev(σᵢⱼ)
    τ = get_τ(sᵢⱼ,r.damage)
    G = r.elasticity.G
    ν = r.elasticity.ν
    
    # stress tensor calculation
    term1 = ( 1 + A1*B1*σ/(2*τ) + B1^2/2) * σoop
    term2 = (3ν/(1+ν) + A1*B1*σ/(2*τ) - A1^2/3 + B1^2/2) * σ
    term3 = A1*B1*τ/3
    return (1/(2G)) * (term1 - term2 + term3)
end
compute_ϵ_oop(r,D,σᵢⱼ)=compute_ϵ_oop(r,compute_A1B1(r,D)...,σᵢⱼ)

function compute_ϵ_oop(r::Rheology,A1,B1,σᵢⱼ::SymmetricTensor)

    σoop = σᵢⱼ[3,3]
    σ = 1/3*tr(σᵢⱼ)
    sᵢⱼ = dev(σᵢⱼ)
    τ = get_τ(sᵢⱼ,r.damage)
    G = r.elasticity.G
    ν = r.elasticity.ν
    
    # stress tensor calculation
    term1 = ( 1 + A1*B1*σ/(2*τ) + B1^2/2) * σoop
    term2 = (3ν/(1+ν) + A1*B1*σ/(2*τ) - A1^2/3 + B1^2/2) * σ
    term3 = A1*B1*τ/3
    return (1/(2G)) * (term1 - term2 + term3)
end

function compute_ϵij(r,A1,B1,σᵢⱼ::SymmetricTensor)

    if all(σᵢⱼ.==0)
        return σᵢⱼ
    end
    σ = 1/3*tr(σᵢⱼ)
    sᵢⱼ = dev(σᵢⱼ)
    τ = get_τ(sᵢⱼ,r.damage)
    G = r.elasticity.G
    ν = r.elasticity.ν
    Id = SymmetricTensor{2,3}(δ)
    # strain invariants 
    term1 = ( 1 + A1*B1*σ/(2*τ) + B1^2/2) * σᵢⱼ
    term2 = (3ν/(1+ν) + A1*B1*σ/(2*τ) - A1^2/3 + B1^2/2) * σ
    term3 = A1*B1*τ/3
    return (1/(2G)) * (term1 + (term3 - term2)*Id)
end

compute_ϵij(r::Rheology,D,σᵢⱼ::SymmetricTensor) = (D==0) ? compute_ϵij(r,0.0,0.0,σᵢⱼ) : compute_ϵij(r,compute_A1B1(r,D)...,σᵢⱼ)
  
function compute_σij(r,A1,B1,Γ,ϵij)
    # TODO make a visco elastic version of this function

    if all(x->x==0,ϵij)
        return SymmetricTensor{2,3}(ϵij)
    end

    G = r.elasticity.G
    ν = r.elasticity.ν
    Id = SymmetricTensor{2,3}(δ)
    # strain invariants
    ϵ = tr(ϵij)
    e = dev(ϵij)
    γ = sqrt(2.0 * e ⊡ e)
    (γ == 0) && (γ = nextfloat(0.0))
    (γ == Inf) && (γ = prevfloat(Inf))

    # stress tensor calculation
    term1 = ( (3*(1-2ν))/(1+ν) + A1^2 - A1*B1*ϵ/γ ) * ϵij
    term2 = (3ν/(1+ν) + B1^2/2 - A1^2/3 + A1*B1*ϵ/(3γ)) * ϵ
    term3 = -A1*B1*γ/2
    return (G/Γ) * (term1 + (term2 + term3)*Id)
end

# CHECKED WITH ALL USED FUNCTIONS
function compute_σij(r,D,ϵij)
    # TODO make a visco elastic version of this function

    if all(x->x==0,ϵij)
        return SymmetricTensor{2,3}(ϵij)
    end
    # unpack
    G = r.elasticity.G
    ν = r.elasticity.ν

    # Create Identity 3*3 tensor
    Id = SymmetricTensor{2,3}(δ)

    # Damage constants
    c1, c2, c3 = compute_c1c2c3(r,D)
    A, B = compute_AB(r,c1,c2,c3)
    A1 = compute_A1(r,A)
    B1 = compute_B1(r,B)
    Γ = compute_Γ(r,A1,B1)

    # strain invariants
    ϵ = tr(ϵij)
    e = dev(ϵij)
    γ = sqrt(2.0 * e ⊡ e)
    (γ == 0) && (γ = nextfloat(0.0))
    (γ == Inf) && (γ = prevfloat(Inf))
    # stress tensor calculation
    term1 = ( (3*(1-2ν))/(1+ν) + A1^2 - A1*B1*ϵ/γ ) * ϵij
    term2 = (3ν/(1+ν) + B1^2/2 - A1^2/3 + A1*B1*ϵ/(3γ)) * ϵ
    term3 = -A1*B1*γ/2
    return (G/Γ) * (term1 + (term2 + term3)*Id)
end

function compute_σ̇ij(r,D,Ḋ,ϵij,ϵ̇ij)

    if all(x->x==0,ϵij)
        all(x->x==0,ϵ̇ij) && (return SymmetricTensor{2,3}(ϵij))
    end
    # unpack
    G = r.elasticity.G
    ν = r.elasticity.ν

    # Create Identity 3*3 tensor
    𝕀 = SymmetricTensor{2,3}(δ)

    # Damage constants
    c1, c2, c3 = compute_c1c2c3(r,D)
    A, B = compute_AB(r,c1,c2,c3)
    A1 = compute_A1(r,A)
    B1 = compute_B1(r,B)
    Γ = compute_Γ(r,A1,B1)
    dc1dD = compute_dc1dD(r,D)
    dc2dD = compute_dc2dD(r,D)
    dc3dD = compute_dc3dD(r,D)
    dA1dD = compute_dA1dD(r,dc1dD,dc2dD,dc3dD,c2,c3)
    dA1dD = compute_dB1dD(r,dc1dD,dc2dD,dc3dD,c2,c3)
    dΓdD = compute_dΓdD(r,A1,B1,dA1dD,dB1dD)

    # strain invariants
    ϵ = tr(ϵij)
    e = dev(ϵij)
    γ = sqrt(2.0 * e ⊡ e)
    (γ == 0) && (γ = nextfloat(0.0))
    (γ == Inf) && (γ = prevfloat(Inf))

    ϵ̇ = tr(ϵ̇ij)
    ė = dev(ϵ̇ij)
    γ̇ = e ⊡ ė / (2*γ)
    
    # stress tensor calculation
    f1 = f₁(ϵ,γ,A1,B1)
    f2 = f₂(ϵ,γ,A1,B1)
    f3 = f₃(A1,B1)
    df1dt = df₁dt(ϵ, γ, A1, B1, ϵ̇, γ̇, Ḋ, dA1dD, dB1dD)
    df2dt = df₂dt(ϵ, γ, A1, B1, ϵ̇, γ̇, Ḋ, dA1dD, dB1dD)
    df3dt = df₃dt(D,Ḋ,dA1dD, dB1dD)

    term1 = df1dt*ϵij + f1*ϵ̇ij + (df2dt*ϵ + f2*ϵ̇ + df3dt*γ + f3*γ̇)*𝕀
    term2 = f1*ϵij + (f2*ϵ + f3*γ)*𝕀
    return (G/Γ) * (term1 - (1/Γ)*dΓdD*Ḋ*term2)
end
f₁(ϵ,γ,A1,B1) = 3*(1-2ν)/(1+ν) + A1^2 - A1*B1*ϵ/γ
f₂(ϵ,γ,A1,B1) = 3*ν/(1+ν) + B1^2/2 - A1^2/3 + A1*B1*ϵ/(3*γ)
f₃(A1,B1) = -A1*B1/2

df₁dt(ϵ, γ, A1, B1, ϵ̇, γ̇, Ḋ, dA1dD, dB1dD) = -(A1*B1/γ)*ϵ̇ + (A1*B1*ϵ/γ^2)*γ̇ + (2*A1*dA1dD - (dA1dD*B1 + A1*dB1dD)*ϵ/γ)*Ḋ
df₂dt(ϵ, γ, A1, B1, ϵ̇, γ̇, Ḋ, dA1dD, dB1dD) = (A1*B1/(3*γ))*ϵ̇ - (A1*B1*ϵ/(3*γ^2))*γ̇ + (B1*dB1dD - 2/3*A1*dA1dD - (dA1dD*B1 + A1*dB1dD)*ϵ/(3γ))*Ḋ
df₃dt(D,Ḋ,dA1dD, dB1dD) = -0.5*(dA1dD*B1 + A1*dB1dD)*Ḋ

function state_system!(du,u,p,t)
    D, ϵ, γ, σ, τ = u
    r, dϵdt, dγdt = p

    # damage constants
    (D < r.damage.D₀) && (println("D = ",D); @warn "D < D0, something went wrong")
    (D == r.damage.D₀) && (D += 1e-9) # insure D > D0 to prevent singularity
    isnan(D) && println("D is NaN")

    c1, c2, c3 = compute_c1c2c3(r,D)
    A, B = compute_AB(r,c1,c2,c3)

    # TODO : Check KI sign to avoid unnecessary calculations

    A1 = compute_A1(r,A)
    B1 = compute_B1(r,B)
    Γ = compute_Γ(r,A1,B1)
    a1 = compute_a1(B1,Γ)
    b1 = compute_b1(A1,B1,Γ)
    b2 = compute_b2(r,A1,Γ)

    # derivatives
    # D
    dDdt = compute_subcrit_damage_rate(r, σ, τ, D)
    #println("D = ", D)
    #println("dDdt = ", dDdt, "\n")

    # c1, c2, c3
    dc1dt = compute_dc1dD(r,D) * dDdt
    dc2dt = compute_dc2dD(r,D) * dDdt
    dc3dt = compute_dc3dD(r,D) * dDdt

    # println("dDdt : ", dDdt)
    # println("dc1dt : ",dc1dt)
    # println("dc2dt : ",dc2dt)
    # println("dc3dt : ",dc3dt)
    # A1, B1
    dA1dt = compute_dA1dt(r,dc1dt,dc2dt,dc3dt,c2,c3)
    dB1dt = compute_dB1dt(r,dc1dt,dc2dt,dc3dt,c2,c3)

    # Γ
    dΓdt = compute_dΓdt(r,A1,B1,dA1dt,dB1dt)

    # a1, b1, b2
    da1dt = compute_da1dt(B1,Γ,dB1dt,dΓdt)
    db1dt = compute_db1dt(A1,B1,Γ,dA1dt,dB1dt,dΓdt)
    db2dt = compute_db2dt(r,A1,Γ,dA1dt,dΓdt)

    # println("Γ : ",Γ)
    # println("A1 : ",A1)
    # println("B1 : ",B1)
    # println("dΓdt : ",dΓdt) # issue here with Nan
    # println("dA1dt : ",dA1dt) # issue here with Nan
    # println("dB1dt : ",dB1dt) # issue here with Nan
    # println("a1,b1,da1dt,db1dt,ϵ,γ,dϵdt,dγdt : ")
    # println(a1)
    # println(b1)
    # println(da1dt) # issue here with Nan
    # println(db1dt) # issue here with Nan
    # println(ϵ)
    # println(γ)
    # println(dϵdt)
    # println(dγdt)
    du[1] = dD = dDdt
    du[2] = dϵ = dϵdt
    du[3] = dγ = dγdt
    du[4] = dσ = compute_dσdt(r,a1,b1,da1dt,db1dt,ϵ,γ,dϵdt,dγdt)
    du[5] = dτ = compute_dτdt(r,b1,b2,db1dt,db2dt,ϵ,γ,dϵdt,dγdt)

    # assertions
    @assert dD >= 0
end

#### TEST ####

# elasticity = Elasticity(E = 70e9,
#                                     ν = 0.3)
# plasticity = DruckerPrager(μ = 0.6,
#                                       C = 1e6)
# damage = BRSDamage( μ = 0.6, # Friction coef
#                     β = 0.1, # Correction factor
#                     K₁c = 1.74e6, # Critical stress intensity factor (Pa.m^(1/2))
#                     a = 1e-3, # Initial flaw size (m)
#                     ψ = atan(0.7), # crack angle to the principal stress (radians)
#                     D₀ = 0.3, #0.1 Initial flaw density
#                     n = 34.0, # Stress corrosion index
#                     l̇₀ = 0.24, # Ref. crack growth rate (m/s)
#                     H = 50e3, # Activation enthalpy (J/mol)
#                     A = 5.71 ) # Preexponential factor (m/s)
# r = Rheology(damage,nothing,elasticity,nothing)
#
# D = range(r.damage.D₀,1.0,length = 50)
# c1 = compute_c1.(Ref(r),D)
# c2 = compute_c2.(Ref(r),D)
# c3 = compute_c3.(Ref(r),D)
# A = compute_A.(Ref(r),c1,c2,c3)
# B = compute_B.(c1,c2,c3)
#
# A1 = compute_A1.(Ref(r),A)
# B1 = compute_B1.(Ref(r),B)
#Γ = compute_Γ.(Ref(r),A1,B1)
#Γ = compute_Γ.(Ref(r),D)
# convex_cond = 1 ./ Γ
#
# p = lineplot(D,convex_cond)
# lineplot!(p,D[1:end-1],A1[1:end-1])
# lineplot!(p,D[1:end-1],Γ[1:end-1])
# show(p)
# pyplot()
#
# σ = -10e6
# τ = 55e6
# D_vec = [0.1:0.001:0.999;]
# KI_vec = compute_KI.(Ref(r), Ref(σ), Ref(τ), D_vec)
# dDdt_vec = compute_subcrit_damage_rate.(Ref(r), KI_vec, D_vec)
# p = plot(D_vec,dDdt_vec)
#     ylims!(0,1)
# display(p)
# p2 = plot(D_vec,KI_vec./r.damage.K₁c)
#     ylims!(0,10)
# display(p2)
