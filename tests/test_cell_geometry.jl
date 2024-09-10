include("../src/PyrochloreGeometry.jl")

import .PyrochloreGeometry as geom
using StaticArrays


# Checking correctness of tetra_idx
#
function verify_all(L::Int, celldelta=[0,0,0])
	dx = 8*L*SVector{3,Int64}(celldelta)
	println("L=$(L), shift=$(dx)")
	lat = geom.PyroFCC(L)
	for (J, t) in enumerate(lat.tetra_sites)
		I = geom.tetra_idx(lat, t+dx)
		@assert I == J """Incorrect index returned.
		Prompt: t=$(t+dx).
		Expected $(J) -> $(t), got $(I) -> $(lat.tetra_sites[I])
		"""
	end
end

for L=1:10
	verify_all(L)
	verify_all(L,[1,0,0])
	verify_all(L,[-1,0,0])
	verify_all(L,[1,1,0])
	verify_all(L,[1,1,1])
	verify_all(L,[1,-1,0])
	verify_all(L,[-1,-1,0])
	verify_all(L,[-1,-1,-1])
	verify_all(L,[1,-1,0])
end
