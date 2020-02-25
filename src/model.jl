using LinearAlgebra

include("/Users/leo/.julia/dev/Rheologies/dev/input_file_traction.jl")
function solve(grid::Grid,variables::Tuple,var_interp_order::Tuple,
               quad_order::Int, quad_type::Symbol,
               bc_dicts::BoundaryConditions,
               rheology::Rheology ;
               filename = nothing)
    dh, dbc, cellvalues_u, cellvalues_p, facevalues_u, K = setup_model(grid::Grid, variables::Tuple,
                                                                    var_interp_order::Tuple,
                                                                    quad_order::Int, quad_type::Symbol,
                                                                    bc_dicts::BoundaryConditions)
    # Assemble stiffness matrix and force vector
    K, f = doassemble(cellvalues_u, cellvalues_p, facevalues_u, K, grid, dh, bc_dicts, rheology);
    apply!(K, f, dbc)
    u = Symmetric(K) \ f;

    # export
    (filename == nothing) && (filename = rheology_summary(rheology))
    J.vtk_grid(filename, dh) do vtkfile
        vtk_point_data(vtkfile, dh, u)
        vtk_point_data(vtkfile, dbc)
    end
    return u
end
