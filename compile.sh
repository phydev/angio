#!/bin/bash
rm -fr mod || echo "Creating mod directory..."
mkdir mod
dir='./src'
ifort -r8 -mkl -O3 -module mod $dir/global.F90 $dir/dderivatives.F90 $dir/misc.F90 $dir/init.F90 $dir/source.F90 $dir/etcdyn.F90 $dir/thinning.F90 $dir/blood_flow.F90 $dir/congraph.F90 $dir/run_angio.F90 $dir/main.F90 -o angio_flow
