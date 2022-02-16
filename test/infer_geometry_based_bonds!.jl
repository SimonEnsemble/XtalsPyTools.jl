module test_infer_geometry_based_bonds!

using Test, Xtals, XtalsPyTools, Graphs

function visual_check(xtal::String)
    c = Crystal(xtal)
    c = replicate(c, (2, 2, 2))
    strip_numbers_from_atom_labels!(c)
    infer_geometry_based_bonds!(c, true) # must
    write_xyz(c, "temp/c.xyz")
    write_bond_information(c, "temp/$(c.name)")
    @info c.name * " see .vtk and .xyz to visually check bonds"
end


if !isdir("temp")
    mkdir("temp")
end

@testset "scipy dep" begin
    scipy = rc[:scipy]
    rc[:scipy] = nothing
    @test_throws ErrorException infer_geometry_based_bonds!(Crystal("IRMOF-1.cif"), true)
    rc[:scipy] = scipy
end

@testset "NiPyC2 Tests" begin
    # NiPyC2 bonding
    c = Crystal("NiPyC2_relax.cif")
    strip_numbers_from_atom_labels!(c)
    c = replicate(c, (4, 4, 4))
    bonding_rules = [BondingRule(:H, :*, 1.2),
                     BondingRule(:Ni, :O, 2.5),
                     BondingRule(:Ni, :N, 2.5),
                     BondingRule(:*, :*, 1.9)]
    infer_bonds!(c, true, bonding_rules=bonding_rules)
    c1 = deepcopy(c)
    conn_comps = connected_components(c.bonds)

    @test length(conn_comps) == 2 # interpenetrated

    c_red = getindex(c, conn_comps[1])
    c_blue = getindex(c, conn_comps[2])

    @test ne(c_red.bonds) + ne(c_blue.bonds) == ne(c.bonds)

    remove_bonds!(c)
    infer_geometry_based_bonds!(c, true)

    @test length(conn_comps) == 2 # interpenetrated

    # debugging outputs (delete me)
    write_bond_information.([c, c1], ["temp/c", "temp/c1"])
    write_xyz(c, "temp/c")

    @test c1.bonds == c.bonds # consistency between two different bonding schemes
end

@testset "FIQCEN Tests" begin
    # FIQCEN bonding
    c = Crystal("FIQCEN_clean.cif")
    strip_numbers_from_atom_labels!(c)
    infer_geometry_based_bonds!(c, true)

    @test length(connected_components(c.bonds)) == 1 # not interpenetrated

    @test c.atoms.species[neighbors(c.bonds, 1)] == [:Cu, :O, :O, :O, :O]

    visual_check("FIQCEN_clean.cif")

    # reduce covalant radius to see Cu-Cu bond disappear
    covalent_radii = rc[:covalent_radii]
    covalent_radii[:Cu] = 1.15
    remove_bonds!(c)

    @test ne(c.bonds) == 0

    infer_geometry_based_bonds!(c, true, covalent_radii=covalent_radii)

    @test c.atoms.species[neighbors(c.bonds, 1)] == [:O, :O, :O, :O]
end

end