#!/bin/bash
echo "Run this script from the tests directory"

mkdir -p testoutput

JFILE="test_gauge_corr.jl"
# JFILE="./test_piflux_indep_uc.jl"
OUTDIR="testoutput"

for b in `echo 0.0 0.2`; do
	#for pt in `echo \\\\Gamma X W L K U`; do
	for pt in `echo \\\\Gamma X X2 W W2 K K2 L U`; do
		echo "Bx=${b}, measuring ${pt} point..."
		julia $JFILE $1 "${pt}" "${b}"
	done
done

echo "check $OUTDIR for output (all curves should be on top of one another)"

