
"Thanks to `koehlerson` and his [Catalyst](https://github.com/koehlerson/Catalyst) for inspiration, and code reuse."
module GmshParser

export extract_mesh

function getNodes(meshString::Array{String})
    nodeStart = findall(meshString .== "\$Nodes")[1] #\$ escape string key $
    nodeEnd = findall(meshString .== "\$EndNodes")[1]
    _, nnodes,min_node_tag,max_node_tag = parse.(Int, split(meshString[nodeStart+1]))
    @assert max_node_tag == nnodes "nodes tags are not contiguous"
    nodes = []
    node = []
    tag = []
    nodeBlock::Int = 0
    n_entity_nodes = 0
    entity_dim = 0
    entity_tag = 0
    tmp = 0
    first = true

    for nodeLine in meshString[nodeStart+2:nodeEnd-1]
        if nodeBlock == 0 && parse(Int, split(nodeLine)[4]) == 0
            continue
        elseif nodeBlock == 0 && parse(Int, split(nodeLine)[4]) != 0
            entity_dim, entity_tag, _, n_entity_nodes = parse.(Int, split(nodeLine))
            #nodeBlock = parse(Int, split(nodeLine)[4])
            nodeBlock = 2 * n_entity_nodes
            tmp = nodeBlock
            if first == false
                #block = push!.(node, tag)
                block = vcat.(tag, node)
                push!(nodes, block)
                tag = []
                node = []
            end
            first = false
        else
            if nodeBlock > tmp // 2
                push!(tag, [entity_dim, entity_tag, parse(Int, nodeLine)])
            else
                push!(node, parse.(Float64, split(nodeLine)))
            end
            nodeBlock -= 1.0
        end
    end
    # The arithmetic above does not account the last node Block
    # hence we need to explicitely append it
    #block = push!.(node, tag)
    block = vcat.(tag, node)
    push!(nodes, block)
    # flatten the array
    return collect(Iterators.flatten(nodes)) # [entDim entTag nodeTag X Y Z]
end

function map_entities_to_physical_names(meshstring)
    linestart = findall(meshstring .== "\$Entities")[1] #\$ escape string key $
    linestop = findall(meshstring .== "\$EndEntities")[1]
    phys_tag2name = map_physical_tags_to_names(meshstring)
    npoints, ncurves, nsurfs, nvols = parse.(Int, split(meshstring[linestart+1]))
    limits = cumsum((npoints, ncurves, nsurfs, nvols))
    entities_tag2name = (Dict{Int,String}(), Dict{Int,String}(), Dict{Int,String}(), Dict{Int,String}())
    dim_id = 1
    for (i,line) in enumerate(meshstring[linestart+2 : linestop-1])
        dim_id = get_dim_id_recursively(i,dim_id,limits)
        if dim_id == 1
            pointline = split(line)
            if length(pointline) > 5
                n_phys_tags = parse(Int,pointline[5])
                for itag in 1:n_phys_tags
                    entity_tag = parse(Int,pointline[1])
                    phys_tag = parse(Int,pointline[5+itag])
                    try
                        push!(entities_tag2name[dim_id], entity_tag => phys_tag2name[phys_tag])
                    catch e
                        if e isa KeyError
                            @info "This KeyError is frequent and due to dicrepancies between defined entities and physical groups names.\nPlease ensure that all entities have a name"
                            throw(KeyError(phys_tag))
                        else
                            throw(e)
                        end
                    end
                end
            end
        elseif dim_id >= 2
            line_vec = split(line)
            if length(line_vec) > 8
                n_phys_tags = parse(Int,line_vec[8])
                for itag in 1:n_phys_tags
                    entity_tag = parse(Int,line_vec[1])
                    phys_tag = parse(Int,line_vec[8+itag])
                    try
                        push!(entities_tag2name[dim_id], entity_tag => phys_tag2name[phys_tag])
                    catch e
                        if e isa KeyError
                            @info "This KeyError is frequent and due to dicrepancies between defined entities and physical groups names.\nPlease ensure that all entities have a name"
                            throw(KeyError(phys_tag))
                        else
                            throw(e)
                        end
                    end
                end
            end
        end
    end
    return entities_tag2name
end

function map_physical_tags_to_names(meshstring)
    linestart = findall(meshstring .== "\$PhysicalNames")[1] #\$ escape string key $
    linestop = findall(meshstring .== "\$EndPhysicalNames")[1]
    phys_tag2name = Dict{Int,String}()
    for line in meshstring[linestart+2 : linestop-1]
        entity_line = split(line)
        phys_tag2name[parse(Int,entity_line[2])] = entity_line[3][2:end-1]
    end
    return phys_tag2name
end

function get_dim_id_recursively(i,dim_id,limits)
    if i > limits[dim_id]
        return get_dim_id_recursively(i,dim_id + 1,limits)
    else
        return dim_id
    end
end

function getElements(meshString::Array{String})
    dim = get_dim(meshString)
    elementStart = findall(meshString .== "\$Elements")[1] #\$ escape string key $
    elementStop = findall(meshString .== "\$EndElements")[1]
    elementTotal = parse(Int, split(meshString[elementStart+1])[2]) #no. of Nodes
    nodeElements = []
    boundaryElements = []
    domainElements = []
    ele = []
    eleBlock::Int = 0
    tmp = 0
    eleDim = 0
    first = true
    eletype = 0
    entity_tag = 1

    for eleLine in meshString[elementStart+2:elementStop-1]
        if eleBlock == 0 && parse(Int, split(eleLine)[4]) == 0 # Entity definition and containing no elements
            continue
        elseif eleBlock == 0 && parse(Int, split(eleLine)[4]) != 0 # entity def and contains elements
            eleBlock = parse(Int, split(eleLine)[4]) # number of element for this entity
            tmp = eleBlock
            if first == false
                if eleDim == dim-1
                    push!(boundaryElements, pushfirst!.(ele,entity_tag))
                elseif eleDim == dim
                    push!(domainElements, ele)
                end
            end
            eleDim, entity_tag, eletype, param = parse.(Int, split(eleLine)) # dimension of the entity
            ele = []
            first = false # TODO move that before loop ends
        else
            push!(ele, parse.(Int, split(eleLine)))  # elem id + nodes global indices
            eleBlock -= 1.0
        end
    end
    if eleDim == dim - 2
        push!(boundaryElements, pushfirst!.(ele,entity_tag))
    elseif eleDim == dim - 1
        push!(boundaryElements, pushfirst!.(ele,entity_tag)) # [ entity_tag  elem id  nodes_gidx... ]
    elseif eleDim == dim
        push!(domainElements, ele) # [ elem id  nodes_idx... ]
    end

    # flatten the array
    nodeElements = collect(Iterators.flatten(nodeElements))
    boundaryElements = collect(Iterators.flatten(boundaryElements))
    domainElements = collect(Iterators.flatten(domainElements))
    return nodeElements, boundaryElements, domainElements
end

# function nodeToDataFrame(nodes)
#     df = DataFrame(tag = [], x = [], y = [], z = [])
#     for node in nodes
#         push!(df, [node[1], node[2], node[3], node[4]])
#     end
#     return df
# end

function get_facesets(meshstring, elements, el_type)
    nodeElements, boundaryElements, domainElements = elements
    n_nodeElements = size(nodeElements,1)
    n_boundaryElements = size(boundaryElements,1)
    n_domainElements, nnodes_domainElements = size(domainElements)
    nnodes_domainElements -= 1
    boundary_faces = Tuple{String,Int,Int}[]
    n_boundaries_entities = length(unique(boundaryElements[:,1]))
    #boundaries_names = get_boundaries_entities_names(meshstring, dim)
    curves_tag2name = map_entities_to_physical_names(meshstring)[2] # get this for curves entities
    for i = 1:n_boundaryElements
        boundary_entity_tag = boundaryElements[i,1]
        boundary_name = curves_tag2name[boundary_entity_tag]
        boundary_nodes = boundaryElements[i,3:end]
        for j = 1:n_domainElements
            element_nodes = domainElements[j,2:end]
            #println("element type is : $(domainElements[j,1]) $el_type")
            if issubset(boundary_nodes, element_nodes)
                faces_nodes_indices = faces_indices(el_type)
                face_id = findall([issubset(element_nodes[face_idxs], boundary_nodes) for face_idxs in faces_nodes_indices])
                @assert length(face_id) == 1 "length(face_id) = $(length(face_id))"
                push!(boundary_faces, (boundary_name, j, face_id[1])) # pushes (physical_entity_tag, elem_id, face_number)
            end # domainElements[j,1] - (n_nodeElements + n_boundaryElements)
        end
    end

    facesets = Dict{String,Set{Tuple{Int,Int}}}()
    for name in unique([face[1] for face in boundary_faces])#unique(two_c[:,1])
        slice = filter(x->x[1]==name, boundary_faces)
        facesets[name] = Set{Tuple{Int,Int}}(getindex.(slice,[2:3]))
    end
    return facesets
end

function get_nodesets(meshstring, nodes) #[entDim entTag nodeTag X Y Z]
    nodes_in_nodalentity = [node for node in nodes if node[1] == 0]#nodes[nodes[:,1] .== 0 ,:]
    nodesets = Dict{String,Set{Int}}()
    points_tag2name = map_entities_to_physical_names(meshstring)[1]
    nnodessets = length(points_tag2name)
    for (entity_tag,entity_name) in pairs(points_tag2name)
        nodes_w_entity_tag = [node for node in nodes_in_nodalentity if node[2] == entity_tag]#[nodes_in_nodalentity[:,2] .== entity_tag,:]
        nodes_tags = Set([node[3] for node in nodes_w_entity_tag])
        push!(nodesets, entity_name => nodes_tags)
    end
    return nodesets
end

function faces_indices(c::String)
    if c in ("Line","QuadraticLine")
        return [1, 2]
    elseif c in ("Triangle","QuadraticTriangle")
        return ([1,2], [2,3], [3,1])
    elseif c in ("Quadrilateral","QuadraticQuadrilateral")
        return ([1,2], [2,3], [3,4], [4,1])
    elseif c in ("Tetrahedron","QuadraticTetrahedron")
        return ([1,2,3], [1,2,4], [2,3,4], [1,4,3])
    elseif c in ("Hexahedron","QuadraticHexahedron")
        return ([1,4,3,2], [1,2,6,5], [2,3,7,6], [3,4,8,7], [1,5,8,4], [5,6,7,8])
    else

    end
end

get_el_type_str(meshstring,higher_order_elem) = _get_el_type_str(get_dim(meshstring),length(higher_order_elem[end,2:end]))

function _get_el_type_str(dim, nnodes)
    if dim == 2
        if nnodes == 3
            el_type = "Triangle" #TRI3
        elseif nnodes == 6
            el_type = "QuadraticTriangle" #TRI3
        elseif nnodes == 4
            el_type = "Quadrilateral" #QUAD4
        elseif nnodes == 9
            el_type = "QuadraticQuadrilateral" #QUAD4
        else
            @error "element type identifier $el_type_number from .msh doesn't have any JuAFEM correspondance ."
        end
    else
        @error "dimension != 2 are not implemented yet"
    end
    return el_type
end

function get_dim(meshstring)
    linestart = findall(meshstring .== "\$Entities")[1] #\$ escape string key $
    ndim_entities = parse.(Int, split(meshstring[linestart+1]))
    domain_dim = findlast(ndim_entities .!= 0) - 1
    return domain_dim
end

function extract_mesh(meshfile)
    io = open(meshfile)
    lines = readlines(io)
    close(io)

    meshnodes = getNodes(lines)

    el_0D, el_1D, el_2D = getElements(lines)

    el_0D_mat = hcat(el_0D...)'
    el_1D_mat = hcat(el_1D...)'
    el_2D_mat = hcat(el_2D...)'
    elements_mat = (el_0D_mat, el_1D_mat, el_2D_mat)

    higher_order_el_mat = elements_mat[findlast(.!isempty.(elements_mat))]
    #println(higher_order_el_mat)
    el_type_str = get_el_type_str(lines,higher_order_el_mat)

    facesets = get_facesets(lines,elements_mat,el_type_str)
    nodesets = get_nodesets(lines, meshnodes)

    return el_type_str, higher_order_el_mat, meshnodes, facesets, nodesets

end

end # module
