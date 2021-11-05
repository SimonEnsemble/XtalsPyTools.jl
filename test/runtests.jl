testfiles = [
    "misc.jl", 
    "infer_geometry_based_bonds!.jl", 
    "primitive_cell.jl"]

using XtalsPyTools, Test

XtalsPyTools.banner()

for testfile ∈ testfiles
    @info "Running test/$testfile"
    include(testfile)
end

@info "Done."