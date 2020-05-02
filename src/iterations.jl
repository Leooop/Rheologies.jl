
function iterate(cellvalues, facevalues, K, grid, dh, bcd, bcn, mp, clock::Clock{Nothing})

    # Doassemble dispatch on the length of cellValues
    rheology = mp[1]
    K, f = doassemble(cellvalues..., facevalues, K, grid, dh, bcn, mp, rheology);

    #Apply Dirichlet boundary conditions
    apply!(K, f, bcd)

    # Solve
    u = Symmetric(K) \ f
end

function iterate(cellvalues, facevalues, K, grid, dh, bcd, bcn, mp, clock::Clock{T}) where {T}

    while clock.current_time <= clock.tspan[2] # TODO
        # Doassemble dispatch on the length of cellValues
        rheology = mp[1]
        K, f = doassemble(cellvalues..., facevalues, K, grid, dh, bcn, mp, rheology);

        #Apply Dirichlet boundary conditions
        apply!(K, f, bcd)

        # Solve
        u = Symmetric(K) \ f;
    end
end
