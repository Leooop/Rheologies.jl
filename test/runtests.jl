using Rheologies ; const R = Rheologies
using JuAFEM ; const J = JuAFEM
using Test
using SafeTestsets

@testset "variables interpolations" begin
    @test R.get_variables_interpolation((:u,), (2,), Cell{2,3,3}) ==
        Lagrange{2,RefTetrahedron,2}()
    @test R.get_variables_interpolation((:u,), (1,), Cell{2,3,3}) ==
        Lagrange{2,RefTetrahedron,1}()
    @test R.get_variables_interpolation((:u,:p), (2,1), Cell{2,3,3}) ==
        (Lagrange{2,RefTetrahedron,2}(), Lagrange{2,RefTetrahedron,1}())
    @test R.get_variables_interpolation((:u,:p), (2,1), Cell{2,4,4}) ==
        (Lagrange{2,RefCube,2}(), Lagrange{2,RefCube,1}())
end

@testset "get quad rules" begin
    @test typeof(R.get_quadrature_rules(2, :legendre, Cell{2,4,4})[1]) ==
        typeof(QuadratureRule{2,RefCube}(2))
    @test typeof(R.get_quadrature_rules(2, :legendre, Cell{2,4,4})[2]) ==
        typeof(QuadratureRule{1,RefCube}(2))
end

@testset "create values" begin
    el_geom = Cell{2,3,3}
    qr, qr_face = R.get_quadrature_rules(3, :legendre, el_geom)
    interp_u1 = R.get_variables_interpolation((:u,), (2,), el_geom)
    interp_u2 = R.get_variables_interpolation((:u,:p), (2,1), el_geom)
    interp_geom = J.default_interpolation(el_geom)
    @test length(R.create_values(qr, qr_face, interp_geom, interp_u1)) == 2
    @test length(R.create_values(qr, qr_face, interp_geom, interp_u2...)) == 3
end


@testset "create material properties" begin
    func_1D(x) = x[1]
    func_2D(x) = x[1] + x[2]
    elasticity1D = IncompressibleElasticity(E = func_1D, ν =0.4999)
    elasticity2D = IncompressibleElasticity(E = func_2D, ν =0.4999)

    r1D = Rheology(nothing,nothing,elasticity1D,nothing)
    r2D = Rheology(nothing,nothing,elasticity2D,nothing)

    grid1D = generate_grid(Line,(2,), Vec(0.0), Vec(10.0))
    grid2D = generate_grid(Quadrilateral,(2, 2), Vec(0.0,0.0), Vec(10.0,10.0))

    mp1D = create_material_properties(grid1D, r1D)
    mp2D = create_material_properties(grid1D, r1D)

    @test r1D(1.0).elasticity.E == 1.0
    @test r2D([1.0,2.0]).elasticity.E == 3.0




end
