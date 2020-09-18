#### OUTPUT WRITER ####

"""
`OutputWriter` contains information about requested model outputs.

# Fields
- `format::Symbol` : exported file format. `:VTK`, `:JLD2` or `:MAT`
- `path::String` : path to the folder
- `outputs::Dict{Symbol,Function}` : dictionary mapping `Symbol`s output names to functions of `r` (rheology) and `s` (state)
- `frequency::F` : save output every `frequency` seconds. To disable frequency sampling use `nothing` as field value
- `interval::I` : save output every `interval` iterations. To disable interval sampling use `nothing` as field value
- `data::Dict{Symbol,Vector{Float64}}` : preallocated output data
- `last_output_time::Float64` : last simulation time where an export occured
- `opened_path::Bool` : `True` if the folder has already been created
"""
mutable struct OutputWriter{TF,F,I}
    path::String
    outputs::Dict{Symbol,Function} # functions of r (rheology) and s (state)
    frequency::F
    interval::I
    data::Dict{Symbol,Vector{Float64}} # preallocation of output data
    last_output_time::Float64
    force_path::Bool
end

# Define aliases for differents file formats
VTKOutputWriter{F,I} = OutputWriter{:VTK,F,I}
JLD2OutputWriter{F,I} = OutputWriter{:JLD2,F,I}
MATOutputWriter{F,I} = OutputWriter{:MAT,F,I}

# Multiple output formats :
struct MixedOutputWriter{OF1,F1,I1,OF2,F2,I2}
    ow1::OutputWriter{OF1,F1,I1}
    ow2::OutputWriter{OF2,F2,I2}
end

"""
    OutputWriter(format, model, path, outputs ; frequency = nothing, interval = nothing)

`OutputWriter` constructor. One of the kwargs `frequency` or `interval` has to be given.

# Positional arguments
- `format::Symbol` : exported file format. `:VTK`, `:JLD2` or `:MAT`
- `model::Model`
- `path::String` : path to the folder
- `outputs::Dict{Symbol,Function}` : dictionary mapping `Symbol`s output names to functions of `r` (rheology) and `s` (state)

# Keyword arguments
- `frequency::F` : save output every `frequency` seconds
- `interval::I` : save output every `interval` iterations
- `force_path::Bool = true` : If false and path already exists, a new path variant is created
"""
function OutputWriter(format, model, path, outputs, frequency, interval, force_path)

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
    TF = typeof(frequency)
    TI = typeof(interval)

    path = create_dir(path, force_path)

    # loading relevant dependencies
    (format == :JLD2) && (@info("loading FileIO package") ; @eval(Rheologies, import FileIO))
    (format == :MAT) && (@info("loading MAT package") ; @eval(Rheologies, import MAT))

    return OutputWriter{format,TF,TI}(path, outputs, frequency, interval, data, model.clock.tspan[1], force_path)
end
OutputWriter(format, model, path, outputs, ::Nothing, ::Nothing, force_path) = @error "Please provide one of the two kwargs `frequency` or `interval`"
OutputWriter(format, model, path, outputs ; frequency = nothing, interval = nothing, force_path = true) = OutputWriter(format, model, path, outputs, frequency, interval, force_path)

VTKOutputWriter(model, path, outputs ; frequency = nothing, interval = nothing, force_path = true) = OutputWriter(:VTK, model, path, outputs, frequency, interval, force_path)
JLD2OutputWriter(model, path, outputs ; frequency = nothing, interval = nothing, force_path = true) = OutputWriter(:JLD2, model, path, outputs, frequency, interval, force_path)
MATOutputWriter(model, path, outputs ; frequency = nothing, interval = nothing, force_path = true) = OutputWriter(:MAT, model, path, outputs, frequency, interval, force_path)


# NOT NEEDED ANYMORE
# function JuAFEM.reinit!(ow::OutputWriter,model)
#     ow.last_output_time = model.clock.tspan[1]
#     ow.opened_path = false
# end
#
# function JuAFEM.reinit!(mow::MixedOutputWriter,model)
#     JuAFEM.reinit!(mow.ow1,model)
#     JuAFEM.reinit!(mow.ow2,model)
# end
### write_output! ###

write_output!(model, u, ::Nothing) = nothing

function write_output!(model, u, ow::OutputWriter{TF,F,I}) where {TF,F,I}
    c = model.clock

    if c.current_time != c.tspan[2]
        if (I <: Real) & (F <: Real)
            c1 = ((model.clock.iter-1)%ow.interval != 0)
            c2 = (model.clock.current_time - ow.last_output_time < ow.frequency)
            !(c1 || c2) && return nothing
        elseif I <: Real
            ((model.clock.iter-1)%ow.interval != 0) && (return nothing)
        elseif F <: Real
            (model.clock.current_time - ow.last_output_time < ow.frequency) && (return nothing)
        end
    end
    # reinitialize data arrays
    for value in values(ow.data)
        fill!(value,0.0)
    end

    # fill data
    nqp = length(model.material_state[1])
    for (el, cell_states) in enumerate(model.material_state)
        r = model.material_properties[el]
        for state in cell_states
            s = state
            for (key,value) in pairs(ow.outputs)
                ow.data[key][el] += value(r,s)/nqp
            end
        end
    end

    # filename
    filename = "iter_"*@sprintf("%04d",model.clock.iter)*"_time_$(model.clock.current_time)"

    # export
    export_sim(filename, model, u, ow)

    # update outputwriter output time
    ow.last_output_time = model.clock.current_time

    return nothing
end

function write_output!(model, u, mixed_ow::MixedOutputWriter)
    write_output!(model, u, mixed_ow.ow1)
    write_output!(model, u, mixed_ow.ow2)
end

# function create_dir!(ow::OutputWriter)
#     path = ow.path
#     if ispath(path)
#         ow.force_path && (return nothing) # we reuse the existing folder
#         if occursin(r"\([0-9]+\)", path[end-4:end])
#             ind1 = findlast('(',path)
#             ind2 = findlast(')',path)
#             num = parse(Int,path[ind1+1:ind2-1])
#             ow.path = path[1:ind1-1]*"("*string(num+1)*")"
#             create_dir!(ow)
#         else
#             ow.path = path*"(1)"
#             create_dir!(ow)
#         end
#     else
#         splitp = splitpath(path)
#         folder = splitp[end]
#         folderdir = joinpath(splitp[1:end-1]...)
#         cd(folderdir)
#         mkdir(folder)
#         ow.opened_path = true
#         return nothing
#     end
# end

function create_dir(path, force_path)
    if ispath(path)
        if force_path
            non_jl_files = filter!(file -> (!occursin(".jl",file) | occursin(".jld",file)),cd(readdir,path))
            if !isempty(non_jl_files)
                @info "directory already contains non .jl files. Do you want to erase them (y/n) ? "
                erase_input(path,non_jl_files)
            end
            return path
        else
            if occursin(r"\([0-9]+\)", path[end-4:end])
                ind1 = findlast('(',path)
                ind2 = findlast(')',path)
                num = parse(Int,path[ind1+1:ind2-1])
                path = path[1:ind1-1]*"("*string(num+1)*")"
                create_dir(path, force_path)
            else
                path = path*"(1)"
                create_dir(path, force_path)
            end
        end
    else
        splitp = splitpath(path)
        folder = splitp[end]
        folderdir = joinpath(splitp[1:end-1]...)
        cd(folderdir)
        mkdir(folder)
        return path
    end
end

function erase_input(path,files)
    key = readline(stdin)
    if (key == "y") | (key == "Y")
        @info "All non .jl files removed from directory"
        rm.(joinpath.(Ref(path),files))
    elseif (key == "n") | (key == "N")
        @info "Existing non .jl files left in place"
    else
        @warn "input key is non recognized, please type \"y\" or \"n\""
        erase_input()
    end
    return nothing
end

function export_sim(filename, model, u, ow::VTKOutputWriter)
    vtk_grid(joinpath(ow.path, filename), model.dofhandler) do vtkfile
        vtk_point_data(vtkfile, model.dofhandler, u) # displacement field
        vtk_point_data(vtkfile, model.dirichlet_bc) # dirichlet boundary conditions
        for (key,value) in pairs(ow.data)
            vtk_cell_data(vtkfile, value, string(key))
        end
    end
    println("VTK saved")
end

function export_sim(filename, model, u, ow::JLD2OutputWriter)
    ow.data[:u] = Vector{eltype(u)}(undef,length(u))
    ow.data[:u] .= u

    # convert keys from symbols to strings
    data_dict = Dict{String,Any}(string(k)=>v  for (k,v) in pairs(ow.data))
    data_dict["t"] = model.clock.current_time # add simulation time
    data_dict["res"] = spatial_resolution(model.grid)
    data_dict["domain_size"] = get_spatial_domain_size(model.grid)

    FileIO.save(joinpath(ow.path, filename*".jld2"), data_dict)
    println("JLD2 saved")
end

function export_sim(filename, model, u, ow::MATOutputWriter)
    ow.data[:u] .= u
    # convert keys from symbols to strings
    data_dict = Dict{String,Any}(string(k)=>v  for (k,v) in pairs(ow.data))
    data_dict["t"] = model.clock.current_time # add simulation time
    data_dict["res"] = [ i for i in spatial_resolution(model.grid)]
    data_dict["domain_size"] = [ i for i in get_spatial_domain_size(model.grid)]
    #save(joinpath(ow.path, filename*".mat"), data_dict)
    println(typeof.(data_dict[key] for key in keys(data_dict)))
    MAT.matwrite(joinpath(ow.path, filename*".mat"), data_dict)
    println("MAT saved")
end

function spatial_resolution(grid::Grid{dim,Cell{2,4,4}}) where {dim}
    lx, ly = get_spatial_domain_size(grid)

    cell = grid.cells[1]
    el_lx = grid.nodes[cell.nodes[2]].x[1] - grid.nodes[cell.nodes[1]].x[1]
    el_ly = grid.nodes[cell.nodes[4]].x[2] - grid.nodes[cell.nodes[1]].x[2]

    return Int(round(lx/el_lx)), Int(round(ly/el_ly))
end

function get_spatial_domain_size(grid)
    x_nodes = [ node.x[1] for node in grid.nodes ]
    y_nodes = [ node.x[2] for node in grid.nodes ]
    lx = maximum(x_nodes) - minimum(x_nodes)
    ly = maximum(y_nodes) - minimum(y_nodes)
    return lx, ly
end
