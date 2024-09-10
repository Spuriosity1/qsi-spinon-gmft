#!/bin/bash
echo "Run this script from the tests directory"


JFILE="test_indep_of_unitcell.jl"
OUTDIR="testoutput"
LOGDIR="tmp"
SIMNAME="indep_uc"

mkdir -p $OUTDIR
mkdir -p $LOGDIR


DATE=`date -I`
for b in `echo 0.0 0.2`; do
	#for pt in `echo \\\\Gamma X W L K U`; do
	for pt in `echo L W \\\\Gamma K X U`; do
		echo "Bx=${b}, measuring ${pt} point..."
		julia $JFILE "${SIMNAME}_${DATE}" "${pt}" "${b}" > "tmp/$pt_$1$2.out" 
	done
done

echo "check $OUTDIR for output (all curves should be on top of one another)"
rm -ri tmp

