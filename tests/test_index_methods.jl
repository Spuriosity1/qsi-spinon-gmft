include("../src/PyrochloreGeometry.jl")

using .PyrochloreGeometry

lattices = [
PyroPrimitive(1,1,1),
PyroPrimitive(1,2,3),
PyroPrimitive(2,3,5),
PyroPrimitive(2,2,2),
PyroPrimitive(5,7,9)
]


function test_tetras(l::PyroPrimitive)
	A = lattice_vectors(l)
	for (J,t) in enumerate(l.tetra_sites)
		@assert tetra_idx(l,t) == J
		@assert tetra_idx(l,t+A*[1,0,0]) == J
		@assert tetra_idx(l,t+A*[-1,0,0]) == J
		@assert tetra_idx(l,t+A*[-1,-1,0]) == J
		@assert tetra_idx(l,t+A*[-1,-1,-1]) == J
		@assert tetra_idx(l,t+A*[1,0,0]) == J
		@assert tetra_idx(l,t+A*[-7,5,-14]) == J
	end
end

failures = 0
for l in lattices
	try			
		test_tetras(l)
	catch e
		failures += 1
		println("Failed test: $(e)")
	end
end

println("Failed $(failures) tests of $(length(lattices))")

