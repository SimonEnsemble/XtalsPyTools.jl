"""
    infer_geometry_based_bonds!(crystal, include_bonds_across_periodic_boundaries::Bool)

Infers bonds by first finding which atoms share a Voronoi face, and then bond the atoms if the distance
 between them is less than the sum of the covalent radius of the two atoms (plus a tolerance).

# Arguments
- `crystal::Crystal`: The crystal structure
- `include_bonds_across_periodic_boundaries::Bool`: Whether to check across the periodic boundaries
- `r::Float`: voronoi radius, Å
- `pad::Float`: amount to pad covalent radii, Å
- `covalent_radii::Dict{Symbol,Float64}`: See [`covalent_radii`](@ref)
"""
function infer_geometry_based_bonds!(crystal::Crystal,
        include_bonds_across_periodic_boundaries::Bool;
        r::Float64=6., pad::Float64=0.,
        covalent_radii::Dict{Symbol,Float64}=rc[:covalent_radii])
    @assert ne(crystal.bonds) == 0 @sprintf("The crystal %s already has bonds. Remove them with the `remove_bonds!` function before inferring new ones.", crystal.name)
    dm = Xtals.distance_matrix(crystal, include_bonds_across_periodic_boundaries)
    for i = 1:crystal.atoms.n
        for j in bonded_atoms(crystal, i, dm, r, pad, covalent_radii)
            add_edge!(crystal.bonds, i, j)
        end
    end
    bond_sanity_check(crystal)
end


"""
    ids_bonded = bonded_atoms(crystal, i, dm; r=6., σ=3.)

Returns the ids of atoms that are bonded to atom `i` by determining bonds using a Voronoi method and covalent radius data (see [`covalent_radii`](@ref))

# Arguments
- `crystal::Crystal`: Crystal structure in which the bonded atoms will be determined
- `i::Int`: Index of the atom we want to determine the bonds of
- `dm::Array{Float64, 2}`: The distance matrix, see [`distance_matrix`](@ref)
- `r::Float64`: The maximum distance used to determine the neighborhood of atom `i`
- `covalent_radii::Dict{Symbol, Float64}`: Cordero parameter dictionary. See [`covalent_radii`](@ref)
- `pad::Float64`: The amount to pad the covalent radius in Å

# Returns
- `ids_bonded::Array{Int, 1}`: A list of indices of atoms bonded to atom `i`
"""
function bonded_atoms(crystal::Crystal, i::Int, dm::Array{Float64, 2},
        r::Float64, pad::Float64, covalent_radii::Dict{Symbol, Float64})
    species_i = crystal.atoms.species[i]
    ids_neighbors, xs, _ = neighborhood(crystal, i, r, dm)
    ids_shared_voro_faces = _shared_voronoi_faces(ids_neighbors, xs)
    ids_bonded = Int[]
    for j in ids_shared_voro_faces
        species_j = crystal.atoms.species[j]
        # sum of covalent radii
        radii_sum = covalent_radii[species_j] + covalent_radii[species_i]
        # margin = σ e.s.d.s, unless that's too small
        max_dist = radii_sum + pad
        if dm[i, j] ≤ max_dist
            push!(ids_bonded, j)
        end
    end
    return ids_bonded
end


"""
    ids_shared_voro_face = _shared_voronoi_faces(ids_neighbors, xs)

Of the neighboring atoms, find those that share a Voronoi face.

# Arguments
- `ids_neighbors::Array{Int, 1}`: indices of atoms within the neighborhood of a specific atom.
- `xs::Array{Array{Float64, 1}, 1}`: array of Cartesian position of the atoms within the neighborhood of a specific atom, relative to the specific atom.

# Returns
- `ids_shared_voro_face::Array{Int, 1}`: indices of atoms that share a Voronoi face with a specific atom
"""
function _shared_voronoi_faces(ids_neighbors::Array{Int,1}, xs::Array{Array{Float64,1},1})
    if isnothing(rc[:scipy])
        error("Python dependency scipy not loaded.")
    else
        scipy = rc[:scipy]
    end
    # first element of xs is the point itself, the origin
    @assert length(ids_neighbors) == (length(xs) - 1)
    voro = scipy.Voronoi(xs)
    rps = voro.ridge_points # connection with atom zero are connection with atom i
    ids_shared_voro_face = Int[] # corresponds to xs, not to atoms of crystal
    for k = 1:size(rps)[1]
        if sort(rps[k, :])[1] == 0 # a shared face with atom i!
            push!(ids_shared_voro_face, sort(rps[k, :])[2])
        end
    end
    # zero based indexing in Scipy accounted for since xs[0] is origin, atom i.
    return ids_neighbors[ids_shared_voro_face]
end


"""
    ids_neighbors, xs, rs = neighborhood(crystal, i, r, dm)

Find and characterize the neighborhood of atom `i` in the crystal `crystal`.
A neighborhood is defined as all atoms within a distance `r` from atom `i`.
The distance matrix `dm` is used to find the distances of all other atoms in the crystal from atom `i`.

# Arguments
- `crystal::Crystal`: crystal structure
- `i::Int`: Index of the atom (in `crystal`) which the neighborhood is to be characterized.
- `r::Float64`: The maximum distance the neighborhood will be characterized.
- `dm::Array{Float64, 2}`: The distance matrix, see [`distance_matrix`](@ref)

# Returns
- `ids_neighbors::Array{Int, 1}`: indices of `crystal.atoms` within the neighborhood of atom `i`.
- `xs::Array{Array{Float64, 1}, 1}`: array of Cartesian positions of the atoms surrounding atom `i`.
    The nearest image convention has been applied to find the nearest periodic image. Also, the coordinates of atom `i`
    have been subtracted off from these coordinates so that atom `i` lies at the origin of this new coordinate system.
    The first vector in `xs` is `[0, 0, 0]` corresponding to atom `i`.
    The choice of type is for the Voronoi decomposition in Scipy.
- `rs::Array{Float64, 1}`: list of distances of the neighboring atoms from atom `i`.
"""
function neighborhood(crystal::Crystal, i::Int, r::Float64, dm::Array{Float64, 2})
    # get indices of atoms within a distance r of atom i
    #  the greater than zero part is to not include itself
    ids_neighbors = findall((dm[:, i] .> 0.0) .& (dm[:, i] .< r))
    # rs is the list of distance of these neighbors from atom i
    rs = [dm[i, id_n] for id_n in ids_neighbors]
    @assert all(rs .< r)
    # xs is a list of Cartesian coords of the neighborhood
    #  coords of atom i are subtracted off
    #  first entry is coords of atom i, the center, the zero vector
    #  remaining entries are neighbors
    # this list is useful to pass to Voronoi for getting Voronoi faces
    #  of the neighborhood.
    xs = [[0.0, 0.0, 0.0]] # this way atom zero is itself
    for j in ids_neighbors
        # subtract off atom i, apply nearest image
        xf = crystal.atoms.coords.xf[:, j] - crystal.atoms.coords.xf[:, i]
        nearest_image!(xf)
        x = crystal.box.f_to_c * xf
        push!(xs, x)
    end
    return ids_neighbors, xs , rs
end