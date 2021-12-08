#!/usr/bin/env bash

gfortran -o angio-flow -fdefault-real-8 src/global.F90 src/dderivatives.F90 src/misc.F90 src/init.F90 src/source.F90 src/etcdyn.F90 src/thinning.F90 src/blood_flow.F90 src/congraph.F90 src/run_angio.F90 src/main.F90 /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

./angio-flow 001