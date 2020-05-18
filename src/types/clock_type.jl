### Time Handler ###
"""
`Clock` type contains information related to time

# Fields
- `tspan`::Tuple{Float64,Float64} : simulation time interval
- `current_time`::Float64 : current simulation time
- `iter`::Int : iteration number
- `Δt_max`::Float64 : maximum timestep
- `Δt`::Float64 : timestep
- `time_vec`::Vector{Float64} : vector containing times for each iteration
- `variable_Δt`::Bool : controls the ability to adapt the timestep
- `Δt_fact_up`::Float64 : timestep increase factor
- `Δt_fact_down`::Float64 : timestep decrease factor
"""
mutable struct Clock
    tspan::Tuple{Float64,Float64}
    current_time::Float64
    iter::Int
    Δt_max::Float64
    Δt::Float64
    time_vec::Vector{Float64}
    variable_Δt::Bool
    Δt_fact_up::Float64 # timestep increase factor
    Δt_fact_down::Float64 # timestep decrease factor
    function Clock(tspan,Δt_max,Δt,variable_Δt,Δt_fact_up,Δt_fact_down)
        current_time = float(tspan[1])
        iter = 0
        time_vec = Float64[tspan[1]]
        (Δt_max < Δt) && (Δt_max = Δt)
        if variable_Δt == false
            Δt_fact_up = 1.0
            Δt_fact_down = 1.0
        end
        new(convert(Tuple{Float64,Float64},tspan),current_time,iter,float(Δt_max),float(Δt),time_vec,variable_Δt,float(Δt_fact_up),float(Δt_fact_down))
    end
end

"""
    Clock(tspan::Tuple, Δt ; <keyword arguments>)

Clock constructor with time span `tspan` and timestep `Δt`.

# Keyword arguments
- `Δt_max` = 1.0 : maximum Δt
- `variable_Δt` = true : is the timestep allowed to change
- `Δt_fact_up` = 1.5 : timestep increase multiplier
- `Δt_fact_down` = 0.5 : timestep decrease multiplier
"""
Clock(tspan, Δt; Δt_max = 1.0, variable_Δt = true, Δt_fact_up = 1.5, Δt_fact_down = 0.5) = Clock(tspan,Δt_max,Δt,variable_Δt,Δt_fact_up,Δt_fact_down)

"""
    Clock(<keyword arguments>)

Clock constructor.

# Keyword arguments
- `tspan::Tuple` (mandatory) : time span
- `Δt` (mandatory) : timestep
- `Δt_max` = 1.0 : maximum Δt
- `variable_Δt` = true : is the timestep allowed to change
- `Δt_fact_up` = 1.5 : timestep increase multiplier
- `Δt_fact_down` = 0.5 : timestep decrease multiplier
"""
Clock(;tspan, Δt_max = 1.0, Δt, variable_Δt = true, Δt_fact_up = 1.5, Δt_fact_down = 0.5) = Clock(tspan,Δt_max,Δt,variable_Δt,Δt_fact_up,Δt_fact_down)


function Base.show(io::IO, ::MIME"text/plain", c::Clock)
        print(io, "⌗  Clock instance\n",
        "├── tspan        : $(c.tspan)\n",
        "├── current_time : $(c.current_time)\n",
        "├── iter         : $(c.iter)\n",
        "├── Δt_max       : $(c.Δt_max)\n",
        "├── Δt           : $(c.Δt)\n",
        "├── time_vec     : $(c.time_vec)\n",
        "├── variable_Δt  : $(c.variable_Δt)\n",
        "├── Δt_fact_up     : $(c.Δt_fact_up)\n",
        "└── Δt_fact_down  : $(c.Δt_fact_down)")
end

function timestep!(c::Clock)
    c.Δt = min(c.Δt, c.Δt_max)
    # handle the case when timestepping overshoots the maximum simulation time
    if (c.current_time + c.Δt) > c.tspan[2]
        c.Δt = c.tspan[2] - c.current_time
        c.current_time = c.tspan[2]
    else
        # don't update first iteration time to start at tspan[1]
        (c.iter != 0) && (c.current_time += c.Δt)
    end
    c.iter += 1
    (c.iter != 0) && push!(c.time_vec,c.current_time)
    return nothing
end

function undo_timestep!(c::Clock)
    c.current_time -= c.Δt
    c.iter -= 1
    pop!(c.time_vec)
    return nothing
end

function JuAFEM.reinit!(c::Clock)
    c.iter = 0
    length(c.time_vec)>=2 && (c.Δt = c.time_vec[2] - c.time_vec[1])
    c.current_time = c.tspan[1]
    c.time_vec = Float64[c.tspan[1]]
    return nothing
end
