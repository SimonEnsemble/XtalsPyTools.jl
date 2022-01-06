testfiles = [
    "misc.jl", 
    "infer_geometry_based_bonds!.jl", 
    "primitive_cell.jl"]

@assert VERSION.major == 1
@assert VERSION.minor ≥ 1

using XtalsPyTools, Test
XtalsPyTools.banner()

for testfile ∈ testfiles
    @info "Running test/$testfile"
    include(testfile)
end

@info "Done."