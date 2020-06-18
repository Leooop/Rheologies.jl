
function RK4(du, u, p, update_du!, Δt::T) where T<:AbstractFloat

    # initialize derivative ponderated sum
    du_sum = zeros(eltype(u),size(u))

    # increment derivative ponderated sum
    update_du!(du,u,p,0.0)
    du_sum .+= du

    update_du!(du,u.+((Δt/2.0).*du),p,0.0)
    du_sum .+= 2.0 .* du

    update_du!(du,u.+(Δt/2.0).*du,p,0.0)
    du_sum .+= 2.0 .* du

    update_du!(du,u.+Δt.*du,p,0.0)
    du_sum .+= du

    # update u with the mean of the ponderated sum over Δt
    return (Δt/6).*du_sum .+ u
end

function timestep_RK4_adaptative(du, u, p, update_du!, Δt::T, e₀::T) where T<: AbstractFloat
    u1 = RK4(du, u, p, update_du!, Δt)
    #println("RK4 u1",u1)
    uint = RK4(du, u, p, update_du!, Δt/2)
    #println("RK4 uint",uint)
    u2 = RK4(du, uint, p, update_du!, Δt/2)
    #println("RK4 u2",u2)

    e = abs(u2[1]-u1[1]) # error calculated on damage only
    if e < e₀
        return u1, Δt, min(Δt*abs(e₀/e),Δt*2)
    else
        timestep_RK4_adaptative(du, u, p, update_du!, float(Δt*abs(e₀/e)^2), e₀)
    end
end

function compute_RK4_adaptative(u0,p,update_du!,tspan,Δt_ini,e₀)
    u_vec = Vector{Vector{Float64}}()
    t_vec = Float64[tspan[1]]
    du = zeros(eltype(u0),size(u0))
    u_next = similar(u0)
    cohesion_loss_flag = false
    push!(u_vec,u0)
    next_Δt = Δt_ini
    while t_vec[end] < tspan[2]
        (t_vec[end]+next_Δt > tspan[2]) && (next_Δt = tspan[2] - t_vec[end]) # prevents time overshoot
        u_next, used_Δt, next_Δt = timestep_RK4_adaptative(du, u_vec[end], p, update_du!,next_Δt,e₀)
        push!(u_vec,u_next)
        push!(t_vec,t_vec[end]+used_Δt)

        # break if D >= 1
        (u_next[1] >= 1) && (cohesion_loss_flag = true ; break)
    end
    return u_vec, t_vec, cohesion_loss_flag
end

function compute_RK4_adaptative_push_version!(u_vec,t_vec,u0,p,update_du!,tspan,Δt_ini,e₀)

    du = zeros(eltype(u0),size(u0))
    u_next = similar(u0)
    cohesion_loss_flag = false
    next_Δt = Δt_ini
    while t_vec[end] < tspan[2]
        (t_vec[end]+next_Δt > tspan[2]) && (next_Δt = tspan[2] - t_vec[end]) # prevents time overshoot
        u_next, used_Δt, next_Δt = timestep_RK4_adaptative(du, u_vec[end], p, update_du!,next_Δt,e₀)
        push!(u_vec,u_next)
        push!(t_vec,t_vec[end]+used_Δt)

        # break if D >= 1
        (u_next[1] >= 1) && (cohesion_loss_flag = true ; break)
    end
    return nothing
end
