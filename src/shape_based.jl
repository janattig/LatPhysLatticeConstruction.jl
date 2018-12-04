################################################################################
#
#   GENERAL INTERFACE FOR CREATING LATTICES INSIDE SHAPES
#   (USING ABSOLUTE COORDINATES)
#
################################################################################



# FALLBACK function
function getLatticeInShape(
            :: Type{L},
            unitcell :: U,
            shape    :: Function,
            origin   :: Integer = 1
        ) :: L where {
            D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},
            U<:AbstractUnitcell{S,B},
            DL,LLS,LLB,SL<:AbstractSite{LLS,DL},BL<:AbstractBond{LLB,0},
            L<:AbstractLattice{SL,BL,U}
        }

    # throw an error as this is not implemented yet
    error("Either function 'getLatticeInShape' is not yet implemented for unitcells of type " * string(U) * ", i.e. N=" *
        string(N) * " / D=" * string(D) * " or you passed the wrong lattice type L = " * string(L))
end



# generate a lattice in 2D in the absolute shape given by the function shape (open boundary conditions)
function getLatticeInShape(
            :: Type{L},
            unitcell :: U,
            shape    :: Function,
            origin   :: Integer = 1
        ) :: L where {
            D,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,2},
            U<:AbstractUnitcell{S,B},
            DL,LLS,LLB,SL<:AbstractSite{LLS,DL},BL<:AbstractBond{LLB,0},
            L<:AbstractLattice{SL,BL,U}
        }


    # Format of checklist item
    #   Tuple{S, NTuple, Int64}
    # denoting
    #   Site - UC Indices (i,j,k,...) - Site Index (in UC)

    # checklist for all sites that are checked if added etc
    checklist = Tuple{S, NTuple{2,Int64}, Int64}[]
    # list for all accepted sites that are checked
    accepted  = Tuple{S, NTuple{2,Int64}, Int64}[]

    # List of all accepted bonds
    bond_list = BL[]


    # push the origin site to the checklist
    push!(
        checklist, (deepcopy(site(unitcell,origin)), (0,0), origin)
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
                # the respective site is not in list, create it to check the position, create a new site object
                site_new = deepcopy(site(unitcell, to(b)))
                # set the new position
                point!(site_new, (uc_to[1] .* a1(unitcell)) .+ (uc_to[2] .* a2(unitcell)) .+ point(site_new))
                # if still inside the shape, add it to the list
                if shape(point(site_new))
                    # format for sites on checklist: ([pos_x, pos_y, pos_z], [uc_i, uc_j, uc_k], index_in_UC)
                    push!(checklist, (
                        site_new,
                        uc_to,
                        to(b)
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

# generate a lattice in 3D in the absolute shape given by the function shape (open boundary conditions)
function getLatticeInShape(
            :: Type{L},
            unitcell :: U,
            shape    :: Function,
            origin   :: Integer = 1
        ) :: L where {
            D,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,3},
            U<:AbstractUnitcell{S,B},
            DL,LLS,LLB,SL<:AbstractSite{LLS,DL},BL<:AbstractBond{LLB,0},
            L<:AbstractLattice{SL,BL,U}
        }


    # Format of checklist item
    #   Tuple{S, NTuple, Int64}
    # denoting
    #   Site - UC Indices (i,j,k,...) - Site Index (in UC)

    # checklist for all sites that are checked if added etc
    checklist = Tuple{S, NTuple{3,Int64}, Int64}[]
    # list for all accepted sites that are checked
    accepted  = Tuple{S, NTuple{3,Int64}, Int64}[]

    # List of all accepted bonds
    bond_list = BL[]


    # push the origin site to the checklist
    push!(
        checklist, (deepcopy(site(unitcell,origin)), (0,0,0), origin)
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
                # the respective site is not in list, create it to check the position, create a new site object
                site_new = deepcopy(site(unitcell, to(b)))
                # set the new position
                point!(site_new, (uc_to[1] .* a1(unitcell)) .+ (uc_to[2] .* a2(unitcell)) .+ (uc_to[3] .* a3(unitcell)) .+ point(site_new))
                # if still inside the shape, add it to the list
                if shape(point(site_new))
                    # format for sites on checklist: ([pos_x, pos_y, pos_z], [uc_i, uc_j, uc_k], index_in_UC)
                    push!(checklist, (
                        site_new,
                        uc_to,
                        to(b)
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
function getLatticeInShape(
            unitcell :: U,
            shape    :: Function,
            origin   :: Integer = 1
        ) :: Lattice{S,Bond{LB,0},U} where {
            D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},
            U<:AbstractUnitcell{S,B}
        }

    # call the general function
    return getLatticeByBondDistance(Lattice{S,Bond{LB,0},U}, unitcell, shape, origin)
end


# export the function
export getLatticeInShape





################################################################################
#
#   SPECIFIC FUNCTIONS FOR CREATING LATTICES INSIDE CERTAIN SHAPES
#   (USING ABSOLUTE COORDINATES & POSSIBILITY OF RELATIVE COORDINATES)
#
#   SUPPORTED SHAPES
#   - box
#   - sphere
#
################################################################################



#-------------
#   BOX
#-------------

# MOST GENERIC INTERFACE (FALLBACK)
function getLatticeInBox(
            :: Type{L},
            unitcell   :: U,
            dimensions :: Vector{<:Real},
            center     :: Vector{<:Real},
            origin     :: Integer = 1
        ) :: L where {
            D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},
            U<:AbstractUnitcell{S,B},
            DL,LLS,LLB,SL<:AbstractSite{LLS,DL},BL<:AbstractBond{LLB,0},
            L<:AbstractLattice{SL,BL,U}
        }

    # throw an error as this is not implemented yet
    error("Either function 'getLatticeInBox' is not yet implemented for unitcells of type " * string(U) * ", i.e. N=" *
        string(N) * " / D=" * string(D) * " or you passed the wrong lattice type L = " * string(L))
end

# MOST GENERIC INTERFACE (2D)
function getLatticeInBox(
            :: Type{L},
            unitcell   :: U,
            dimensions :: Vector{<:Real},
            center     :: Vector{<:Real},
            origin     :: Integer = 1
        ) :: L where {
            N,LS,LB,S<:AbstractSite{LS,2},B<:AbstractBond{LB,N},
            U<:AbstractUnitcell{S,B},
            LLS,LLB,SL<:AbstractSite{LLS,2},BL<:AbstractBond{LLB,0},
            L<:AbstractLattice{SL,BL,U}
        }

    # create the shape function
    shape_function(p) = (center[1] - dimensions[1]/2 <= p[1] <= center[1] + dimensions[1]/2) &&
                        (center[2] - dimensions[2]/2 <= p[2] <= center[2] + dimensions[2]/2)

    # call the interface with the shape function
    return getLatticeInShape(L, unitcell, shape_function, origin)
end

# MOST GENERIC INTERFACE (3D)
function getLatticeInBox(
            :: Type{L},
            unitcell   :: U,
            dimensions :: Vector{<:Real},
            center     :: Vector{<:Real},
            origin     :: Integer = 1
        ) :: L where {
            N,LS,LB,S<:AbstractSite{LS,3},B<:AbstractBond{LB,N},
            U<:AbstractUnitcell{S,B},
            LLS,LLB,SL<:AbstractSite{LLS,3},BL<:AbstractBond{LLB,0},
            L<:AbstractLattice{SL,BL,U}
        }

    # create the shape function
    shape_function(p) = (center[1] - dimensions[1]/2 <= p[1] <= center[1] + dimensions[1]/2) &&
                        (center[2] - dimensions[2]/2 <= p[2] <= center[2] + dimensions[2]/2) &&
                        (center[3] - dimensions[3]/2 <= p[3] <= center[3] + dimensions[3]/2)

    # call the interface with the shape function
    return getLatticeInShape(L, unitcell, shape_function, origin)
end


# export the function
export getLatticeInBox


# Wrapper function (2D)
function getLatticeInBox(
            unitcell   :: U,
            dimensions :: Vector{<:Real},
            center     :: Vector{<:Real} = [0.0, 0.0],
            origin     :: Integer = 1
        ) :: Lattice{S,Bond{LB,0},U} where {
            N,LS,LB,S<:AbstractSite{LS,2},B<:AbstractBond{LB,N},
            U<:AbstractUnitcell{S,B}
        }

    # call the interface
    return getLatticeInBox(Lattice{S,Bond{LB,0},U}, unitcell, dimensions, center, origin)
end

# Wrapper function (3D)
function getLatticeInBox(
            unitcell   :: U,
            dimensions :: Vector{<:Real},
            center     :: Vector{<:Real} = [0.0, 0.0, 0.0],
            origin     :: Integer = 1
        ) :: Lattice{S,Bond{LB,0},U} where {
            N,LS,LB,S<:AbstractSite{LS,3},B<:AbstractBond{LB,N},
            U<:AbstractUnitcell{S,B}
        }

    # call the interface
    return getLatticeInBox(Lattice{S,Bond{LB,0},U}, unitcell, dimensions, center, origin)
end
