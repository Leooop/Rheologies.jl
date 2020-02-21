using Rheologies ; const R = Rheologies
using JuAFEM
using Test

@testset "interface.jl" begin
    @test R.get_variables_interpolation((:u,), (2,), Cell{2,3,3}) ==
        Lagrange{2,RefTetrahedron,2}()
    @test R.get_variables_interpolation((:u,), (1,), Cell{2,3,3}) ==
        Lagrange{2,RefTetrahedron,1}()
    @test R.get_variables_interpolation((:u,:p), (2,1), Cell{2,3,3}) ==
        (Lagrange{2,RefTetrahedron,2}(), Lagrange{2,RefTetrahedron,1}())
    @test R.get_variables_interpolation((:u,:p), (2,1), Cell{2,4,4}) ==
        (Lagrange{2,RefCube,2}(), Lagrange{2,RefCube,1}())

    @test typeof(R.get_quadrature_rule(2, :legendre, Cell{2,4,4})) ==
        typeof(QuadratureRule{2,RefCube}(2))
    @test typeof(R.get_quadrature_rule(3, :legendre, Cell{2,3,3})) ==
        typeof(QuadratureRule{2,RefTetrahedron}(3))

end
