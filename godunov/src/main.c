/* ======================================== 
 *  
 * An exact Riemann solver.
 *
 * Usage: ./riemann paramfile.txt ic.dat
 *  
 * ======================================== */



#include <stdlib.h>
#include <stdio.h>

#include "gas.h"
#include "godunov.h"
#include "io.h"
#include "params.h"
#if RIEMANN==EXACT
#include "riemann-exact.h"
#endif

double gamma = 1.4;   /* Ratio of specific heats; Adiabatic exponent... */

double* x = 0;
state* w_new = 0;
state* w_old = 0;
state* w_intercell = 0; /* W[i] corresponds to W[i-1/2] */

double t = 0;
double dt = 0;
double dx = 1;
double vmax = 0;

params pars;



/* ====================================== */
int main(int argc, char* argv[]){
/* ====================================== */

  /* ------------------ */
  /* Initial setup      */
  /* ------------------ */

  init_params();

  read_cmdlineargs(argc, argv);
  read_paramfile();
  if (pars.verbose) print_params();
  state left, right;
  read_ic(&left, &right);



  /* allocate memory for state arrays */
  /* leave extra space for boundaries */
  /* initialize values                */

  x = malloc((pars.nx+4)*sizeof(double));
  dx = 2./(pars.nx);
  for (int i=0; i<pars.nx+4; i++){
    x[i] = i*dx - (1+2*dx);
  }

  w_old = malloc((pars.nx+4)*sizeof(state));
  w_new = malloc((pars.nx+4)*sizeof(state));
  w_intercell = malloc((pars.nx+4)*sizeof(state));

  for (int i=0; i<pars.nx/2+2; i++){
    w_old[i].rho = left.rho;
    w_old[i].u   = left.u;
    w_old[i].p   = left.p;
  }
  for (int i=pars.nx/2+2; i<pars.nx+4; i++){
    w_old[i].rho = right.rho;
    w_old[i].u   = right.u;
    w_old[i].p   = right.p;
  }


  for (int i=0; i<pars.nx+4; i++){
    printf("%10.4lf ", w_old[i].rho);
  }
  printf("\n");
  for (int i=0; i<pars.nx+4; i++){
    printf("%10.4lf ", x[i]);
  }
  printf("\n");
  /* for (int i=0; i<pars.nx+4; i++){ */
  /*   printf("%10.4lf ", w_old[i].p); */
  /* } */
  /* printf("\n"); */

  int step = 0;
  while(t < pars.tmax){

    step += 1;
    if (step == pars.nsteps) break;

    compute_intercell_states();

    /* compute new dt */
    dt = pars.ccfl * dx/vmax;

for (int i=0; i<pars.nx+4; i++){
  printf("%10.4lf ", w_intercell[i].rho);
}
printf("\n");

    t += dt;
    if (pars.verbose) printf("Finished step %d at t = %12.6lf\n", step, t);

  }



/*   check_ic(&left, &right); */
  /*  */
  /*  */
  /* if (!vacuum) { */
  /*   compute_star_state(&left, &right, &starL, &starR); */
  /*   print_results(&left, &right, &starL, &starR); */
  /* } */
  /*  */
  /* if (p.nsteps > 0) { */
  /*  */
  /*    */
  /*   for (int istep=0; istep<p.nsteps+1; istep++){ */
  /*     if (p.verbose) printf("Computing and writing output for step=%d\n", istep); */
  /*     double t = ((double) istep)/((double) p.nsteps)*p.tmax; */
  /*     if (vacuum) { */
  /*       compute_solution_vacuum(t, X, RHO, U, P, &pars, &left, &right); */
  /*     } */
  /*     else{ */
  /*       compute_solution(t, X, RHO, U, P, &pars, &left, &right, &starL, &starR); */
  /*     } */
  /*     write_output(istep, t, X, RHO, U, P, &pars); */
  /*   } */
  /*  */
  /*  */
  /* } */
/*  */

  return(0);

}




