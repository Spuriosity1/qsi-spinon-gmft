module PyrochloreGeometry

export PyroPrimitive 
export tetra_idx, tetra_IDX
export spin_sl, spin_idx
export get_hexagons, get_dual_fcc_locations
export A_sites, B_sites
export lattice_vectors, reciprocal_basis

using StaticArrays

const Vec3 = SVector{3,Int64};
const Vec3_F64 = SVector{3, Float64};


const MVec3 = MVector{3,Int64};
const MVec3_F64 = MVector{3, Float64};



const pyro = map(x->SVector{3,Int}(x), [
     [ 1,  1,  1],
     [ 1, -1, -1],
     [-1,  1, -1],
     [-1, -1,  1]
     ])

# r = [pyro0]*0.125  pyro1]*0.125 
#                               pyro2]*0.125  pyro3]*0.125]

const diamond = map(x->SVector{3,Int}(x),[
    [0, 0, 0], [2, 2, 2]
])

const fcc_Dy = map(x->SVector{3,Int}(x),[
    [0, 0, 0],
    [0, 4, 4],
    [4, 0, 4],
    [4, 4, 0]
    ])

const fcc_Ti = map(x->SVector{3,Int}(x),[
    [4, 4, 4],
    [4, 0, 0],
    [0, 4, 0],
    [0, 0, 4]
    ])


const plaqt = map(y->map(x->SVector{3,Int}(x),y),[
    [
        [ 0, -2,  2],
        [ 2, -2,  0],
        [ 2,  0, -2],
        [ 0,  2, -2],
        [-2,  2,  0],
        [-2,  0,  2]],
    [
        [ 0,  2, -2],
        [ 2,  2,  0],
        [ 2,  0,  2],
        [ 0, -2,  2],
        [-2, -2,  0],
        [-2,  0, -2]
    ],
    [
        [ 0, -2, -2],
        [-2, -2,  0],
        [-2,  0,  2],
        [ 0,  2,  2],
        [ 2,  2,  0],
        [ 2,  0, -2]],
    [
        [ 0,  2,  2],
        [-2,  2,  0],
        [-2,  0, -2],
        [ 0, -2, -2],
        [ 2, -2,  0],
        [ 2,  0,  2]
    ]
    ])


const S2 = sqrt(2)
const S3 = sqrt(3)
const S6 = sqrt(6)

# Local axes [sl][1,2,3]
const axis = map(x->SMatrix{3,3,Float64}(reduce(hcat,x)),[
    [[ 1, 1,-2]./S6, [-1, 1, 0]./S2, [ 1, 1, 1]./S3],
    [[ 1,-1, 2]./S6, [-1,-1, 0]./S2, [ 1,-1,-1]./S3],
    [[-1, 1, 2]./S6, [ 1, 1, 0]./S2, [-1, 1,-1]./S3],
    [[-1,-1,-2]./S6, [ 1,-1, 0]./S2, [-1,-1, 1]./S3]
    ])

"""
Represents a pyrochlore lattice within cubic FCC superlattice
@member L the number of cubic cells along one dimension (e.g. L=2 has 2^3*16 =128 sites) 
@member tetra_sites the tetrahedron locations
@member spin_sites the spin locations
"""
struct PyroFCC
    L::Int
    tetra_sites::Vector{Vec3}
    # A_sites::Vector{Vec3}
    spin_sites::Vector{Vec3}
    
    function PyroFCC(L::Int)
        as = Vector{Vec3}()
        ss = Vector{Vec3}()
        
        for ix in 0:L-1, iy in 0:L-1, iz in 0:L-1
            for fcc=1:4
                rA = fcc_Dy[fcc] + 8*[ix, iy, iz]
                push!(as, rA)
                for psl=1:4
                    push!(ss, rA + pyro[psl])
                end
            end
        end
        new(L, vcat(as, map(x->x+[2, 2, 2], as)), ss)
    end
end

const primitive_basis = @SMatrix [ 0 4 4; 4 0 4; 4 4 0 ]
const primitive_recip_basis = SMatrix{3,3}( 2π/8 * [ -1 1 1; 1 -1 1; 1 1 -1 ])

"""
Represents a pyrochlore lattice with generic basis
@member L the number of primitive cells along the three primitive axes
@member tetra_sites the tetrahedron locations (the first half of threse are 'A' sites
@member spin_sites the spin locations
"""
struct PyroPrimitive
    L::SVector{3, Int}
    tetra_sites::Vector{Vec3} # the first half of these are the 'A' sites
    spin_sites::Vector{Vec3}

    function PyroPrimitive(L1, L2, L3)
        L = @SVector [ L1, L2, L3 ]
        
        as = Vector{Vec3}()
        ss = Vector{Vec3}()
        
        for ix in 0:L[1]-1, iy in 0:L[2]-1, iz in 0:L[3]-1
            rA = primitive_basis * [ix, iy, iz]
            # for fcc=1:4
            #     rA = fcc_Dy[fcc] + 8*[ix, iy, iz]
            push!(as, rA)
            for spin_sl=1:4
                push!(ss, rA + pyro[spin_sl])
            end
        end
        new(L, vcat(as, map(x->x+[2, 2, 2], as)), ss)
    end
end

const PyroGeneric = Union{PyroFCC, PyroPrimitive}

function A_sites(lat::PyroGeneric)
    L = div( length(lat.tetra_sites), 2)
    return lat.tetra_sites[1:L]
end


function B_sites(lat::PyroGeneric)
    L = div( length(lat.tetra_sites), 2)
    return lat.tetra_sites[L+1:end]
end





"""
tetra_idx(lattice::PyroFCC, tetra_pos_::Vec3)

Returns the index of site tetra_pos_ in lattice.tetra_sites.
"""
function tetra_idx(lattice::PyroFCC, tetra_pos::SVector{3,Int64})
    
    diamond_sl =  (tetra_pos[1] & 0x2) >> 1
    #  !!! 0-based index
    
    # if this bit is set, there is a "+2" or "+6" somewhere. 
    # therefore it is SL 2
    # tetra_pos_ -= diamond[diamond_sl]

    L = lattice.L

    # arcane bit math nonsense
    # Relies on the specific form of fcc_Dy. In 0-based, fcc_Dy is
    # 0=> [0,0,0]
    # 1=> [0,4,4]
    # 2=> [4,0,4]
    # 3=> [4,4,0]
    # The first two indices can then be read as bit representations of the index.
    # Since the sum is of the form N*8 + fcc_Dy + [either 2 or 0], 
    # the 0x4 bits of the first two dimensions uniquely identify the sublattice.
    fcc_sl = ( (tetra_pos[2]&0x4) >> 2) | ( (tetra_pos[1]&0x4) >> 1)

    # read off the cell
    #I = mod.(div.(tetra_pos, 8), L)
    I = div.(mod.(tetra_pos, 8L), 8)
    
    
    idx =  (
        (
            (
                diamond_sl*L + I[1]
            )*L + I[2]
        )*L + I[3]
    )*4 + fcc_sl + 1

    @assert all(mod.(lattice.tetra_sites[idx] - tetra_pos, 8L) .==0)
    return idx
end

# old implementation
#=
function tetra_idx(lattice::PyroFCC, tetra_pos_::Vec3)
    L = lattice.L
   
    tmp = MVector{3, Int64}(tetra_pos_)
    tmp2 = MVector{3, Int64}(tetra_pos_)

    #decide the plquette sublattice

    diamond_sl = 2
    # decide the diamond sublattice
    if all(mod.(tmp, 4) .== 0)
        diamond_sl = 1
    else
        tmp -= diamond[2]
    end
    # @assert all(mod.(tmp, 4) .== 0)
    
    fcc_sl = -1
    for fcc in 1:4
        tmp2 = tmp-fcc_Dy[fcc]
        if all( mod.(tmp2, 8) .== 0)
            tmp=tmp2
            fcc_sl = fcc
            break
        end
    end

    I = mod.(div.(tmp,8),L)

    idx =  (diamond_sl-1) * L^3*4 + I[1]*L^2*4 + I[2]*L*4 + I[3]*4 + fcc_sl
    # make sure we did a good job
    # @assert lattice.tetra_sites[idx] == tetra_pos_
   
    return idx
end
=#




@inline function spin_sl(lattice::PyroFCC, spin_pos_::Vec3)
    
    # decide the FCC sl
    
    L = lattice.L
   
    tmp = MVector{3, Int64}(spin_pos_)
    tmp2 = MVector{3, Int64}(spin_pos_)

    #decide the spin sublattice
    spin_sl = -1
    for ssl in 1:4
        tmp2 = tmp - pyro[ssl]
        
        if all(mod.(tmp2, 4) .== 0)
            tmp = tmp2
            spin_sl = ssl
            break
        end
    end
    
    return spin_sl
end

function spin_idx(lattice::PyroFCC, spin_pos_::Vec3)
    
    # decide the FCC sl
    
    L = lattice.L
   
    tmp = MVector{3, Int64}(spin_pos_)
    tmp2 = MVector{3, Int64}(spin_pos_)

    #decide the spin sublattice
    spin_sl = -1
    for ssl in 1:4
        tmp2 = tmp - pyro[ssl]
        
        if all(mod.(tmp2, 4) .== 0)
            tmp = tmp2
            spin_sl = ssl
            break
        end
    end
    
    @assert spin_sl != -1
    fcc_sl = -1
    for fcc in 1:4
        tmp2 = tmp - fcc_Dy[fcc]
        if all(mod.(tmp2, 8) .== 0)
            tmp = tmp2
            fcc_sl = fcc
            break
        end
    end
    @assert fcc_sl != -1
    
    I = mod.(div.(tmp,8),lattice.L)

    idx =  I[1]*L^2*16 + I[2]*L*16 + I[3]*16 + (fcc_sl-1)*4 + spin_sl
    return idx
end


function get_dual_A_fcc_locations(lattice::PyroFCC)
    locs = []
    for ix in 0:lattice.L-1, iy in 0:lattice.L-1, iz in 0:lattice.L-1
        cell = 8 .*[ix,iy,iz]
        for fcc in fcc_Ti
            push!(locs, cell + fcc)
        end
    end
    return locs
end


function get_dual_fcc_locations(lattice::PyroFCC)
    locs = get_dual_A_fcc_locations(lattice)
    append!(locs, map(x->x+[2,2,2], locs))
    return locs
end





const u = @SMatrix [ 0 1 0; 1 0 0; 1 1 -1 ]
const v = @SMatrix [ 1 0 -1; 0 1 -1; 0 0 1 ]
const d = @SVector [ 4, 4, 8 ]

function tetra_IDX(lattice::PyroPrimitive, tetra_pos_::SVector{3,Int64})

    diamond_sl =  (tetra_pos_[1] & 0x2) >> 1
    #  !!! 0-based index
    
    # if this bit is set, there is a "+2" or "+6" somewhere. 
    # therefore it is SL 2
    @assert mod.(u * tetra_pos_, d) == (@SVector [0,0,0]) ||  mod.(u * tetra_pos_, d) == (@SVector [2,2,2])
    return mod.(v * fld.(u * tetra_pos_, d), lattice.L), diamond_sl
end


using LinearAlgebra

function is_Bravais(lat::PyroPrimitive, x,tol=1e-10)
    return norm(exp.(1im .* reciprocal_basis(lat)'*x).- 1) < tol
end

function tetra_idx(lattice::PyroPrimitive, tetra_pos_::SVector{3,Int64})
    I, diamond_sl = tetra_IDX(lattice, tetra_pos_)
    res= ((diamond_sl*lattice.L[1]+I[1])*lattice.L[2]+I[2])*lattice.L[3] + I[3] + 1
    @assert begin
    is_Bravais(lattice, lattice.tetra_sites[res] - tetra_pos_)     
    end "Incorrect lattice resolution: Got $(tetra_pos_), assigned $(lattice.tetra_sites[res])"
    return res
end

function spin_sl(spin_pos_::SVector{3, Int})
    # figure out the spin sublattice
   
    tmp = MVector{3, Int64}(spin_pos_)

    #decide the spin sublattice
    for ssl in 1:4
        tmp = spin_pos_ - pyro[ssl]
        
        if all(mod.(tmp, 4) .== 0)
            return ssl
        end
    end
    
end

function spin_idx(lattice::PyroPrimitive, spin_pos_::SVector{3, Int})
    
    tmp = MVector{3, Int64}(spin_pos_)

    #decide the spin sublattice
    spin_sl = -1
    for ssl in 1:4
        tmp = spin_pos_ - pyro[ssl]
        
        if all(mod.(tmp, 4) .== 0)
            spin_sl = ssl
            break
        end
    end
    
    I, _ = tetra_IDX(lattice, tmp)
    return ((I[1]*lattice.L[2]+I[2])*lattice.L[3] + I[3])*4 + spin_sl
end


function get_dual_A_fcc_locations(lattice::PyroPrimitive)
    locs = []
    for ix in 0:lattice.L[1]-1, iy in 0:lattice.L[2]-1, iz in 0:lattice.L[3]-1
        rA = primitive_basis * [ix, iy, iz]
        push!(locs, rA + [4,4,4])
    end
        
    return locs
end


function get_hexagons(lattice::PyroGeneric)
    #=
    Returns a list of n_spin hexagons, in the format of 6-membered lists
    [[J1, mu1] [J1, nu1] [J2, mu3] .... [J3,nu3]]
    where the J's are 'A' tetrahedron indices and the mu's are pyrochlore sublattices.
    =#
    hexa_sites = [] # index by hexa_sites[dual_fcc_idx][sl] -> 6-member arr
    for (I, r) in enumerate(get_dual_A_fcc_locations(lattice))
        row = []
        for mu in 1:4
            r_plaq = r + pyro[mu]
            spin_sites = [r_plaq+y for y in plaqt[mu] ]

            push!(row, [])
            for (j, s) in enumerate(spin_sites)
                nu = spin_sl(s)
                fcc = s - pyro[nu]
                J = tetra_idx(lattice, fcc)
                push!(last(row), (J, nu))
            end
        end     
        push!(hexa_sites, row)
    end
    return hexa_sites
end

"""
lattice_vectors(lat::PyroPrimitive)
Returns as 3x3 matrix of the basis vectors, stored as columns
"""
function lattice_vectors(lat::PyroPrimitive)
    return lat.L' .* primitive_basis 
end

"""
reciprocal_basis(lat::PyroPrimitive)
Returns as 3x3 matrix of the reciprocal basis vectors, stored as columns
satisfies reciprocal_basis' * lattice_vectors = 2π eye(3)
"""
function reciprocal_basis(lat::PyroPrimitive)
    return inv(lattice_vectors(lat)').*2π
end

function wrap_BZ(lat::PyroPrimitive, Q::SVector{3,Float64})
    B = reciprocal_basis(lat)
    Binv = lattice_vectors(lat)'/(2π)
	return B *( mod.(Binv * Q , 1) )
end

const high_symmetry_points = Dict(
    "\\Gamma"=> [0.,0.,0.],
    "X"=> [1.,0.,0.],
    "W"=> [1.,0.5,0.],
    "K"=> [0.75,0.75,0.],
    "L"=> [0.5,0.5,0.5],
    "U"=> [1.0, 0.25,0.25]
)
end
