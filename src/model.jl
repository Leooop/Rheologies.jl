using LinearAlgebra

#include("/Users/leo/.julia/dev/Rheologies/dev/input_file_traction.jl")
function run_simulation(grid::Grid,variables::Variables,
               quad_order::Int, quad_type::Symbol,
               bc_dicts::BoundaryConditions,
               rheology::Rheology,
               clock::Clock;
               filename = nothing)

    # setup types
    dh, dbc, cellvalues, facevalues, K = setup_model(grid::Grid,
                                                    variables::Variables,
                                                    quad_order::Int, quad_type::Symbol,
                                                    bc_dicts::BoundaryConditions)
    # Assemble stiffness matrix and force vector
    u = assemble_and_solve(cellvalues, facevalues, K, grid, dh, dbc, bc_dicts.neumann, rheology, clock)


    # export
    (filename == nothing) && (filename = rheology_summary(rheology))
    J.vtk_grid(filename, dh) do vtkfile
        vtk_point_data(vtkfile, dh, u)
        vtk_point_data(vtkfile, dbc)
    end
    return u
end
