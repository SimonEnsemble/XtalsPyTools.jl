module test_primitive_cell

using Test, Xtals, XtalsPyTools

@testset "pymatgen dep" begin
    pymatgen = rc[:pymatgen]
    rc[:pymatgen] = nothing
    @test_throws ErrorException primitive_cell(Crystal("IRMOF-1.cif"))
    rc[:pymatgen] = pymatgen
end

@testset "primitive_cell" begin
    xtal = Crystal("str_m1_o3_o20_pcu_sym.22.cif")
    prim = primitive_cell(xtal)
    @test prim.atoms.n == 996
    @test isnothing(assert_P1_symmetry(prim))
end

end
