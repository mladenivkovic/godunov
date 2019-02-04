#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "params.h"


extern double gamma;
extern params pars;


/* ========================================== */
void init_params(){
/* ========================================== */
  /* Initialize parameters to default values  */
  /*------------------------------------------*/

  pars.verbose = 0;
  for (int i = 0; i<MAX_FNAME_SIZE; i++){
    pars.paramfilename[i] = 0;
    pars.datafilename[i] = 0;
  }
  pars.nsteps = 0;
  pars.tmax = 1;
  pars.nx = 100;
  pars.ccfl = 0.9;
  pars.twostate_ic = 1;

}




/* ========================================== */
void print_params(){
/* ========================================== */
  /* Print out current parameters             */
  /*------------------------------------------*/

  printf("==================================================\n");
  printf("Starting calculation. Parameters are:\n");

  printf("Verbose?                 ");
  if (pars.verbose) {
    printf("True\n");
  } else {
    printf("False\n");
  }
  printf("tmax:                    %g\n", pars.tmax);
  printf("nsteps:                  %d\n", pars.nsteps);
  printf("nx:                      %d\n", pars.nx);
  printf("IC has only two states?  ");
  if (pars.twostate_ic) {
    printf("True\n");
  } else {
    printf("False\n");
  }
  printf("gamma:                   %g\n", gamma);

  printf("==================================================\n");

}
