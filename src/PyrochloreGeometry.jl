module PyrochloreGeometry

using StaticArrays

const Vec3 = SVector{3,Int64};
const Vec3_F64 = SVector{3, Float64};

pyro = map(x->SVector{3,Int}(x), [
     [ 1,  1,  1],
     [ 1, -1, -1],
     [-1,  1, -1],
     [-1, -1,  1]
     ])

# r = [pyro0]*0.125  pyro1]*0.125 
#                               pyro2]*0.125  pyro3]*0.125]

diamond = map(x->SVector{3,Int}(x),[
    [0, 0, 0], [2, 2, 2]
])

fcc_Dy = map(x->SVector{3,Int}(x),[
    [0, 0, 0],
    [0, 4, 4],
    [4, 0, 4],
    [4, 4, 0]
    ])

fcc_Ti = map(x->SVector{3,Int}(x),[
    [4, 4, 4],
    [4, 0, 0],
    [0, 4, 0],
    [0, 0, 4]
    ])


plaqt = map(y->map(x->SVector{3,Int}(x),y),[
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


S2 = sqrt(2)
S3 = sqrt(3)
S6 = sqrt(6)

# Local axes [sl][1,2,3]
axis = map(x->SMatrix{3,3,Float64}(reduce(hcat,x)),[
    [[ 1, 1,-2]./S6, [-1, 1, 0]./S2, [ 1, 1, 1]./S3],
    [[ 1,-1, 2]./S6, [-1,-1, 0]./S2, [ 1,-1,-1]./S3],
    [[-1, 1, 2]./S6, [ 1, 1, 0]./S2, [-1, 1,-1]./S3],
    [[-1,-1,-2]./S6, [ 1,-1, 0]./S2, [-1,-1, 1]./S3]
    ])


struct PyroFCC
    L::Int
    tetra_sites::Vector{Vec3}
    A_sites::Vector{Vec3}
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
        new(L, vcat(as, map(x->x+[2, 2, 2], as)), as, ss)
    end
end

function tetra_idx(lattice::PyroFCC, tetra_pos_::Vec3)
    
    # decide the FCC sl
    
    L = lattice.L
   
    tmp = MVector{3, Int64}(tetra_pos_)

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
   
    return idx
end

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


function get_hexagons(lattice::PyroFCC)
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
                nu = spin_sl(lattice, s)
                fcc = s - pyro[nu]
                J = tetra_idx(lattice, fcc)
                push!(last(row), (J, nu))
            end
        end     
        push!(hexa_sites, row)
    end
    return hexa_sites
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
