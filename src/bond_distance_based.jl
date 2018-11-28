#=
# generate a lattice in 3D by bond distance (open boundary conditions)
function getLatticeByBondDistance(
            unitcell        :: U,
            bonddistance    :: Integer,
            origin          :: Integer = 1
        )


    # load the data from the unit cell
    uc_basis            = unitcell.basis
    uc_connections      = unitcell.connections
    uc_lattice_vectors  = unitcell.lattice_vectors

    # arrays for new positions and connections
    positions   = []
    connections = Array{Any, 1}[]

    # checklist for all sites that are checked if added etc
    checklist   = []
    # format for sites on checklist: ([pos_x, pos_y, pos_z], [uc_i, uc_j, uc_k], index_in_UC, bd_current)

    # push the origin site to the checklist
    push!(
        checklist, (uc_basis[origin], [0,0,0], origin, 0)
    )

    # iterate while the checklist is not empty
    while size(checklist,1) > 0

        # get the item that was on the checklist for longest
        item_to_handle = shift!(checklist)

        # check if the item is already in the positions list
        found = false
        for p in positions
            if item_to_handle[2] == p[2] && item_to_handle[3] == p[3]
                # if yes, continue
                found = true
                break
            end
        end
        if found
            continue
        end

        # if not, push it into
        push!(positions, item_to_handle)
        # for all connections
        index_from = size(positions, 1)

        # insert all connections to sites that are already inside the positions list and all other sites into the checklist
        for c in uc_connections
            # check if correct connections
            if c[1] != item_to_handle[3]
                continue
            end
            # search for the other element of the connection
            i_to = item_to_handle[2][1] + c[4][1]
            j_to = item_to_handle[2][2] + c[4][2]
            k_to = item_to_handle[2][3] + c[4][3]
            a_to = c[2]
            index_to = -1
            for (index,item_handled) in enumerate(positions)
                if item_handled[2] == [i_to, j_to, k_to] && item_handled[3] == a_to
                    # make a new connection
                    index_to = index
                    # break the loop
                    break
                end
            end
            # determine whether the element is already inside the list
            if index_to > 0
                # element is at index index_to an can be linked
                connection_new = [index_from; index_to; c[3]; (0, 0, 0)]
                # check if the connetion is already added
                if connection_new in connections
                    continue
                end
                # register as connection
                push!(connections, connection_new)
            else
                # element not in list yet, maybe should be added
                if item_to_handle[4] < bonddistance
                    # format for sites on checklist: ([pos_x, pos_y, pos_z], [uc_i, uc_j, uc_k], index_in_UC, bd_current)
                    push!(checklist, (
                        (i_to * uc_lattice_vectors[1]) .+ (j_to * uc_lattice_vectors[2]) .+ (k_to * uc_lattice_vectors[3]) .+ uc_basis[c[2]],
                        [i_to, j_to, k_to],
                        c[2],
                        item_to_handle[4]+1
                    ))
                end
            end
        end
    end

    # change the format of positions
    positions_TMP = positions

    # erase positions
    positions = Array{Float64, 1}[]
    positions_indices = Int64[]

    # insert the real positions
    for p in positions_TMP
        push!(positions, p[1])
        push!(positions_indices, p[3])
    end

    # insert missing connections (if (i to j) is present, insert (j to i))
    for c in connections
        c_proposed = [c[2]; c[1]; c[3]; (0, 0, 0)]
        if !(c_proposed in connections)
            push!(connections, c_proposed)
        else
            #println("connection already there")
        end
    end

    # save everything to a Lattice object
    lattice = Lattice(
            unitcell,
            Int64[],
            Array{Float64, 1}[],
            positions,
            positions_indices,
            connections,
            filename
        )
    # save the lattice object
    if save
        saveLattice(lattice)
    end

    # return the lattice object
    return lattice
end
=#
