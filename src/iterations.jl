
function assemble_and_solve(cellvalues, facevalues, K, grid, dh, bcd, bcn, rheology, clock::Clock{Nothing}) where {E}

    # Doassemble dispatch on the length of cellValues
    K, f = doassemble(cellvalues..., facevalues, K, grid, dh, bcn, rheology);

    #Apply Dirichlet boundary conditions
    apply!(K, f, bcd)

    # Solve
    u = Symmetric(K) \ f;
end
