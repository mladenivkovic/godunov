#!/bin/bash

make
if [[ $? -ne 0 ]]; then
    exit 1
fi
echo "############################################################################"


genparamfile() {
    # generate parameter file.
    # $1 = nsteps
    # $2 = tmax

    f=paramfile.txt
    echo "// parameter file for godunov program" > $f
    echo ""             >> $f
    echo "verbose = 1"  >> $f
    echo "nx = 100"    >> $f
    echo "gamma = 1.4"  >> $f
    echo "ccfl = 0.1"   >> $f
    echo "nsteps = $1"  >> $f
    echo "tmax = $2"    >> $f
    echo "foutput = $3" >> $f
}

genparamfile 0 1 20
rm -r sod_test
./godunov paramfile.txt ../ic/sod_test.dat
../plot_godunov_solution.py sod_test

genparamfile 0 1 20
rm -r sod_test_modified
./godunov paramfile.txt ../ic/sod_test_modified.dat
../plot_godunov_solution.py sod_test_modified

genparamfile 0 0.5 20
rm -r 123problem
./godunov paramfile.txt ../ic/123problem.dat
../plot_godunov_solution.py 123problem

genparamfile 0 0.025 20
rm -r left_blast_wave
./godunov paramfile.txt ../ic/left_blast_wave.dat
../plot_godunov_solution.py left_blast_wave

genparamfile 0 0.1 20
rm -r right_blast_wave
./godunov paramfile.txt ../ic/right_blast_wave.dat
../plot_godunov_solution.py right_blast_wave

genparamfile 0 0.1 20
rm -r two_shocks
./godunov paramfile.txt ../ic/two_shocks.dat
../plot_godunov_solution.py two_shocks

# genparamfile 0 0.2 100
# rm -r left_vacuum
# ./godunov paramfile.txt ../ic/left_vacuum.dat
# ../plot_godunov_solution.py left_vacuum

# genparamfile 0 0.2 100
# rm -r right_vacuum
# ./godunov paramfile.txt ../ic/right_vacuum.dat
# ../plot_godunov_solution.py right_vacuum

# genparamfile 0 0.2 100
# rm -r vacuum_generating
# ./godunov paramfile.txt ../ic/vacuum_generating.dat
# ../plot_godunov_solution.py vacuum_generating

