#!/bin/bash
echo "Run this script from the tests directory"

mkdir -p testoutput

JFILE="test_indep_of_unitcell.jl"
# JFILE="./test_piflux_indep_uc.jl"
OUTDIR="testoutput"

for b in `echo 0.0 0.2`; do
	#for pt in `echo \\\\Gamma X W L K U`; do
	for pt in `echo L W \\\\Gamma K X U`; do
		echo "Bx=${b}, measuring ${pt} point..."
		julia $JFILE $1 "${pt}" "${b}" $2 > "tmp/$pt_$1$2.out" 
	done
done

echo "check $OUTDIR for output (all curves should be on top of one another)"


