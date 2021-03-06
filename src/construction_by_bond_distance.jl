# FALLBACK function
function getLatticeByBondDistance(
            :: Type{L},
            unitcell        :: U,
            bonddistance    :: Integer,
            origin          :: Integer = 1
        ) :: L where {
            D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},
            U<:AbstractUnitcell{S,B},
            DL,LLS,LLB,SL<:AbstractSite{LLS,DL},BL<:AbstractBond{LLB,0},
            L<:AbstractLattice{SL,BL,U}
        }

    # throw an error as this is not implemented yet
    error("Either function 'getLatticeByBondDistance' is not yet implemented for unitcells of type " * string(U) * ", i.e. N=" *
        string(N) * " / D=" * string(D) * " or you passed the wrong lattice type L = " * string(L))
end



# generate a lattice in 2D by bond distance (open boundary conditions)
function getLatticeByBondDistance(
            :: Type{L},
            unitcell        :: U,
            bonddistance    :: Integer,
            origin          :: Integer = 1
        ) :: L where {
            D,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,2},
            U<:AbstractUnitcell{S,B},
            DL,LLS,LLB,SL<:AbstractSite{LLS,DL},BL<:AbstractBond{LLB,0},
            L<:AbstractLattice{SL,BL,U}
        }


    # Format of checklist item
    #   Tuple{S, NTuple, Int64, Int64}
    # denoting
    #   Site - UC Indices (i,j,k,...) - Site Index (in UC) - Bonddistance to origin

    # checklist for all sites that are checked if added etc
    checklist = Tuple{S, NTuple{2,Int64}, Int64, Int64}[]
    # list for all accepted sites that are checked
    accepted  = Tuple{S, NTuple{2,Int64}, Int64, Int64}[]

    # List of all accepted bonds
    bond_list = BL[]


    # push the origin site to the checklist
    push!(
        checklist, (deepcopy(site(unitcell,origin)), (0,0), origin, 0)
    )

    # iterate while the checklist is not empty
    while size(checklist,1) > 0

        # get the item that was on the checklist for longest
        item_to_handle = popfirst!(checklist)


        # check if the item is already in the accepted list
        found = false
        for acc in accepted
            # check if a site with index acc[3] in unitcell acc[2] is already found
            if item_to_handle[2] == acc[2] && item_to_handle[3] == acc[3]
                found = true
                break
            end
        end
        # if the element is found, continue with the next checklist item
        if found
            continue
        end



        # Item is NOT in the accepted list, push it into the list
        push!(accepted, item_to_handle)
        # set the current item index to the length of the list
        index_from = length(accepted)

        # search through all connections of the unitcell to find the ones relevant persuing
        for b in bonds(unitcell)

            # check if the bond actually describes a bond emerging from the site
            if from(b) != item_to_handle[3]
                continue
            end

            # build up the element that the bond connects to

            # find out in which copy of unitcell the site is located
            uc_to = (
                item_to_handle[2][1] + wrap(b)[1],
                item_to_handle[2][2] + wrap(b)[2]
            )
            # find out which site index it is
            alpha_to = to(b)


            # find out if the item already exists within the accepted list
            index_to = -1
            for (index,item_handled) in enumerate(accepted)
                # check if the item is the particular item
                if item_handled[2] == uc_to && item_handled[3] == alpha_to
                    # save the list index
                    index_to = index
                    # break the loop
                    break
                end
            end
            # determine whether the element is already inside the list (index was set to > 0)
            if index_to > 0
                # create a new bonds
                bond_new_1 = newBond(
                    BL,
                    index_from,
                    index_to,
                    label(b),
                    NTuple{0,Int64}()
                )
                bond_new_2 = newBond(
                    BL,
                    index_to,
                    index_from,
                    label(b),
                    NTuple{0,Int64}()
                )
                # check if the bonds are already added
                if !(bond_new_1 in bond_list)
                    push!(bond_list, bond_new_1)
                end
                if !(bond_new_2 in bond_list)
                    push!(bond_list, bond_new_2)
                end
            else
                # the respective site is not in list, maybe add it to the checklist?
                if item_to_handle[4] < bonddistance
                    # it can be added, create a new site object
                    site_new = deepcopy(site(unitcell, to(b)))
                    # set the new position
                    point!(site_new, (uc_to[1] .* a1(unitcell)) .+ (uc_to[2] .* a2(unitcell)) .+ point(site_new))
                    # format for sites on checklist: ([pos_x, pos_y, pos_z], [uc_i, uc_j, uc_k], index_in_UC, bd_current)
                    push!(checklist, (
                        site_new,
                        uc_to,
                        to(b),
                        item_to_handle[4]+1
                    ))
                end
            end
        end
    end

    # insert missing connections (if (i to j) is present, insert (j to i))
    # TODO

    # save everything to a Lattice object and return it
    return newLattice(
            L,
        	Vector{Float64}[],
        	SL[newSite(SL, point(s[1]), label(s[1])) for s in accepted],
        	bond_list,
            unitcell
        )
end

# generate a lattice in 3D by bond distance (open boundary conditions)
function getLatticeByBondDistance(
            :: Type{L},
            unitcell        :: U,
            bonddistance    :: Integer,
            origin          :: Integer = 1
        ) :: L where {
            D,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,3},
            U<:AbstractUnitcell{S,B},
            DL,LLS,LLB,SL<:AbstractSite{LLS,DL},BL<:AbstractBond{LLB,0},
            L<:AbstractLattice{SL,BL,U}
        }


    # Format of checklist item
    #   Tuple{S, NTuple, Int64, Int64}
    # denoting
    #   Site - UC Indices (i,j,k,...) - Site Index (in UC) - Bonddistance to origin

    # checklist for all sites that are checked if added etc
    checklist = Tuple{S, NTuple{3,Int64}, Int64, Int64}[]
    # list for all accepted sites that are checked
    accepted  = Tuple{S, NTuple{3,Int64}, Int64, Int64}[]

    # List of all accepted bonds
    bond_list = BL[]


    # push the origin site to the checklist
    push!(
        checklist, (deepcopy(site(unitcell,origin)), (0,0,0), origin, 0)
    )

    # iterate while the checklist is not empty
    while size(checklist,1) > 0

        # get the item that was on the checklist for longest
        item_to_handle = popfirst!(checklist)


        # check if the item is already in the accepted list
        found = false
        for acc in accepted
            # check if a site with index acc[3] in unitcell acc[2] is already found
            if item_to_handle[2] == acc[2] && item_to_handle[3] == acc[3]
                found = true
                break
            end
        end
        # if the element is found, continue with the next checklist item
        if found
            continue
        end



        # Item is NOT in the accepted list, push it into the list
        push!(accepted, item_to_handle)
        # set the current item index to the length of the list
        index_from = length(accepted)

        # search through all connections of the unitcell to find the ones relevant persuing
        for b in bonds(unitcell)

            # check if the bond actually describes a bond emerging from the site
            if from(b) != item_to_handle[3]
                continue
            end

            # build up the element that the bond connects to

            # find out in which copy of unitcell the site is located
            uc_to = (
                item_to_handle[2][1] + wrap(b)[1],
                item_to_handle[2][2] + wrap(b)[2],
                item_to_handle[2][3] + wrap(b)[3]
            )
            # find out which site index it is
            alpha_to = to(b)


            # find out if the item already exists within the accepted list
            index_to = -1
            for (index,item_handled) in enumerate(accepted)
                # check if the item is the particular item
                if item_handled[2] == uc_to && item_handled[3] == alpha_to
                    # save the list index
                    index_to = index
                    # break the loop
                    break
                end
            end
            # determine whether the element is already inside the list (index was set to > 0)
            if index_to > 0
                # create a new bonds
                bond_new_1 = newBond(
                    BL,
                    index_from,
                    index_to,
                    label(b),
                    NTuple{0,Int64}()
                )
                bond_new_2 = newBond(
                    BL,
                    index_to,
                    index_from,
                    label(b),
                    NTuple{0,Int64}()
                )
                # check if the bonds are already added
                if !(bond_new_1 in bond_list)
                    push!(bond_list, bond_new_1)
                end
                if !(bond_new_2 in bond_list)
                    push!(bond_list, bond_new_2)
                end
            else
                # the respective site is not in list, maybe add it to the checklist?
                if item_to_handle[4] < bonddistance
                    # it can be added, create a new site object
                    site_new = deepcopy(site(unitcell, to(b)))
                    # set the new position
                    point!(site_new, (uc_to[1] .* a1(unitcell)) .+ (uc_to[2] .* a2(unitcell)) .+ (uc_to[3] .* a3(unitcell)) .+ point(site_new))
                    # format for sites on checklist: ([pos_x, pos_y, pos_z], [uc_i, uc_j, uc_k], index_in_UC, bd_current)
                    push!(checklist, (
                        site_new,
                        uc_to,
                        to(b),
                        item_to_handle[4]+1
                    ))
                end
            end
        end
    end

    # insert missing connections (if (i to j) is present, insert (j to i))
    # TODO

    # save everything to a Lattice object and return it
    return newLattice(
            L,
        	Vector{Float64}[],
        	SL[newSite(SL, point(s[1]), label(s[1])) for s in accepted],
        	bond_list,
            unitcell
        )
end




# wrapper
function getLatticeByBondDistance(
            unitcell        :: U,
            bonddistance    :: Integer,
            origin          :: Integer = 1
        ) :: Lattice{S,Bond{LB,0},U} where {
            D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},
            U<:AbstractUnitcell{S,B}
        }

    # call the general function
    return getLatticeByBondDistance(Lattice{S,Bond{LB,0},U}, unitcell, bonddistance, origin)
end


# export the function
export getLatticeByBondDistance
