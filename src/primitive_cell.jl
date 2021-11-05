"""
    ```prim = primitive_cell(xtal)```

Returns the minimal representation (primitive unit cell) of a crystal structure.
"""
function primitive_cell(xtal::Crystal)
    # check for pymatgen.io.cif dependency
    if isnothing(rc[:pymatgen])
        error("Python dependency pymatgen not loaded.")
    else
        pymatgen = rc[:pymatgen]
    end
    tempfile = joinpath(pwd(), ".temp_$(uuid1()).cif")
    # copy out xtal and convert it to primitive cell
    write_cif(xtal, tempfile)
    pymatgen.CifParser(tempfile).get_structures()[1].get_primitive_structure().to(filename=tempfile)
    # load the result
    primitive = Crystal(tempfile)
    # clean up
    rm(tempfile)
    name = split(xtal.name, ".cif")[1] * "_primitive_cell.cif"
    return Crystal(name, primitive.box, primitive.atoms, primitive.charges)
end
