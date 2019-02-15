#!/bin/bash


#===================================
genparamfile() {
#===================================
    # generate parameter file.
    # $1 = nsteps
    # $2 = tmax
    # $3 = output frequency: after how many steps to write

    f=paramfile.txt
    echo "// parameter file for godunov program" > $f
    echo ""             >> $f
    echo "verbose = 1"  >> $f
    echo "nx = 100"     >> $f
    echo "gamma = 1.4"  >> $f
    echo "ccfl = 0.1"   >> $f
    echo "nsteps = $1"  >> $f
    echo "tmax = $2"    >> $f
    echo "foutput = $3" >> $f
}



#---------------------------------------------
# for SOLVER in HLL; do
# for SOLVER in TSRS; do
for SOLVER in EXACT; do
# for SOLVER in EXACT TRRS; do
# for SOLVER in EXACT TRRS TSRS; do
#---------------------------------------------

    sed -i "s/^RIEMANN=.*/RIEMANN=${SOLVER}/" Makefile

    make clean
    make
    if [[ $? -ne 0 ]]; then
        exit 1
    fi
    echo "############################################################################"



    # genparamfile 1 1 20
    # rm -r $SOLVER/sod_test
    # ./godunov paramfile.txt ../ic/sod_test.dat
    # gdb --args ./godunov paramfile.txt ../ic/sod_test.dat
    # valgrind ./godunov paramfile.txt ../ic/sod_test.dat
    # ../plot_godunov_solution.py $SOLVER/sod_test

    genparamfile 0 1 0
    rm -r $SOLVER/sod_test
    ./godunov paramfile.txt ../ic/sod_test.dat
    ../plot_godunov_solution.py $SOLVER/sod_test

    genparamfile 0 0.1 0
    rm -r $SOLVER/sod_test_modified
    ./godunov paramfile.txt ../ic/sod_test_modified.dat
    ../plot_godunov_solution.py $SOLVER/sod_test_modified

    genparamfile 0 0.2 0
    # genparamfile 30 0.2 5
    rm -r $SOLVER/123problem
    ./godunov paramfile.txt ../ic/123problem.dat
    ../plot_godunov_solution.py $SOLVER/123problem

    genparamfile 0 0.025 0
    rm -r $SOLVER/left_blast_wave
    ./godunov paramfile.txt ../ic/left_blast_wave.dat
    ../plot_godunov_solution.py $SOLVER/left_blast_wave

    genparamfile 0 0.1 0
    rm -r $SOLVER/right_blast_wave
    ./godunov paramfile.txt ../ic/right_blast_wave.dat
    ../plot_godunov_solution.py $SOLVER/right_blast_wave

    genparamfile 0 0.1 0
    rm -r $SOLVER/two_shocks
    ./godunov paramfile.txt ../ic/two_shocks.dat
    ../plot_godunov_solution.py $SOLVER/two_shocks
    #
    # # # genparamfile 0 0.2 0
    # # rm -r $SOLVER/left_vacuum
    # # ./godunov paramfile.txt ../ic/left_vacuum.dat
    # # ../plot_godunov_solution.py $SOLVER/left_vacuum
    #
    # # genparamfile 0 0.2 0
    # # rm -r $SOLVER/right_vacuum
    # # ./godunov paramfile.txt ../ic/right_vacuum.dat
    # # ../plot_godunov_solution.py $SOLVER/right_vacuum
    #
    # # genparamfile 0 0.2 0
    # # rm -r $SOLVER/vacuum_generating
    # # ./godunov paramfile.txt ../ic/vacuum_generating.dat
    # # ../plot_godunov_solution.py $SOLVER/vacuum_generating
    #
    #
done;
