module XtalsPyTools

using Xtals, PyCall, FIGlet, Printf, Graphs, UUIDs
using PrecompileSignatures: @precompile_signatures

# - lists python dependencies
# - checks for python dependencies at pre-compile, warns on failure
# - provides function for python dependency instantiation
include("pydeps.jl")

include("infer_geometry_based_bonds!.jl")
include("primitive_cell.jl")
include("misc.jl")

function __init__()
    # load Python dependencies
    init_pydeps()
    # check Python dependencies
    for pydep ∈ PYDEPS
        check_pydep(pydep)
    end
end

export infer_geometry_based_bonds!, primitive_cell, rc

@precompile_signatures(XtalsPyTools)

end
