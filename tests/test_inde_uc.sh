#!/bin/bash
echo "Run this script from the tests directory"

mkdir -p testoutput

JFILE="test_indep_of_unitcell.jl"
OUTDIR="testoutput"

for b in `echo 0.0 0.02 0.04 0.06`; do
	for pt in `echo \\\\Gamma X W L K U`; do
		echo "Bx=${b}, measuring ${pt} point..."
		julia $JFILE production "${pt}" "${b}"
	done
done

echo "check $OUTDIR for output (all curves should be on top of one another)"


