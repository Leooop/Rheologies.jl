using LinearAlgebra

#include("/Users/leo/.julia/dev/Rheologies/dev/input_file_traction.jl")
function run_simulation(grid::Grid,variables::Variables,
               quad_order::Int, quad_type::Symbol,
               bc_dicts::BoundaryConditions,
               rheology::Rheology,
               clock::Clock;
               filename = nothing)

    # setup types
    dh, bcd, cellvalues, facevalues, mp, K = setup_model(grid::Grid,
                                                    variables::Variables,
                                                    quad_order::Int, quad_type::Symbol,
                                                    bc_dicts::BoundaryConditions,
                                                    rheology::Rheology)
    println(typeof(mp[1]))
    # Assemble stiffness matrix and force vector
    u = iterate(cellvalues, facevalues, K, grid, dh, bcd, bc_dicts.neumann, mp, clock)


    # export
    (filename == nothing) && (filename = rheology_summary(rheology))
    J.vtk_grid(filename, dh) do vtkfile
        vtk_point_data(vtkfile, dh, u)
        vtk_point_data(vtkfile, bcd)
    end
    return u
end
