#### OUTPUT WRITER ####
abstract type AbstractOutputWriter end

mutable struct VTKOutputWriter{F,I} <: AbstractOutputWriter
    path::String
    outputs::Dict{Symbol,Function} # functions of r (rheology) and s (state)
    frequency::F
    interval::I
    data::Dict{Symbol,Vector{Float64}} # preallocation of output data
    last_output_time::Float64
    opened_path::Bool
end

function VTKOutputWriter(model, path, outputs, frequency, interval)
    ncells = getncells(model.grid)
    data = Dict{Symbol,Vector{Float64}}()
    for (key,value) in pairs(outputs)
        data[key] = Vector{Float64}(undef,ncells)
    end
    if frequency isa Nothing
        @assert interval isa Real
    else
        @assert interval isa Nothing
    end
    return VTKOutputWriter(path, outputs, frequency, interval, data, model.clock.tspan[1], false)
end

VTKOutputWriter(model, path, outputs ; frequency = nothing, interval = nothing) = VTKOutputWriter(model, path, outputs, frequency, interval)

function JuAFEM.reinit!(ow::AbstractOutputWriter,model)
    ow.last_output_time = model.clock.tspan[1]
    ow.opened_path = false
end

write_output!(model,::Nothing) = nothing

function write_output!(model, u, ow::VTKOutputWriter{F,I}) where {F,I}
    if I <: Real
        ((model.clock.iter-1)%ow.interval != 0) && (return nothing)
    elseif F <: Real
        (model.clock.current_time - ow.last_output_time < ow.frequency) && (return nothing)
    end
    # reinitialize data arrays
    for value in values(ow.data)
        fill!(value,0.0)
    end

    # fill data
    for (el, cell_states) in enumerate(model.material_state)
        r = model.material_properties[el]
        for state in cell_states
            s = state
            for (key,value) in pairs(ow.outputs)
                ow.data[key][el] += value(r,s)/4
            end
        end
    end

    # export VTK
    ow.opened_path || create_dir!(ow)
    iter_path = joinpath(ow.path, "iter_$(model.clock.iter)_time$(model.clock.current_time)")
    vtk_grid(iter_path, model.dofhandler) do vtkfile
        vtk_point_data(vtkfile, model.dofhandler, u) # displacement field
        vtk_point_data(vtkfile, model.dirichlet_bc) # dirichlet boundary conditions
        for (key,value) in pairs(ow.data)
            vtk_cell_data(vtkfile, value, string(key))
        end
    end
    println("VTK saved")

    ow.last_output_time = model.clock.current_time

    return nothing
end

function create_dir!(ow::AbstractOutputWriter)
    path = ow.path
    if ispath(path)
        if occursin(r"\([0-9]+\)", path[end-4:end])
            ind1 = findlast('(',path)
            ind2 = findlast(')',path)
            num = parse(Int,path[ind1+1:ind2-1])
            ow.path = path[1:ind1-1]*"("*string(num+1)*")"
            create_dir!(ow)
        else
            ow.path = path*"(1)"
            create_dir!(ow)
        end
    else
        splitp = splitpath(path)
        folder = splitp[end]
        folderdir = joinpath(splitp[1:end-1]...)
        cd(folderdir)
        mkdir(folder)
        ow.opened_path = true
        return nothing
    end
end
