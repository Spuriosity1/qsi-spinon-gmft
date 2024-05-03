include("../src/BZMath.jl")
using .BZMath

high_symmetry_points = Dict(
    "\\Gamma"=> [0.,0.,0.],
    "X"=> [1.,0.,0.],
    "W"=> [1.,0.5,0.],
    "K"=> [0.75,0.75,0.],
    "L"=> [0.5,0.5,0.5],
    "U"=> [1.0, 0.25,0.25]
	)


path = generate_path(
	high_symmetry_points, 
	split("\\Gamma X W K \\Gamma L U W"),
	points_per_unit=30,
	K_units=2π/8
	)


