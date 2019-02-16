/* ======================================== 
 *  
 *  A 1D hydrodynamics solver using 
 *  Godunov's method and various Riemann
 *  solvers.
 *  
 * ======================================== */


#include <stdlib.h>
#include <stdio.h>

#include "gas.h"
#include "godunov.h"
#include "io.h"
#include "params.h"

#ifdef RIEMANN_EXACT
#include "riemann-exact.h"
#elif defined RIEMANN_TRRS
#include "riemann-trrs.h"
#elif defined RIEMANN_TSRS
#include "riemann-tsrs.h"
#elif defined RIEMANN_HLL
#include "riemann-hll.h"
#elif defined RIEMANN_HLLC
#include "riemann-hllc.h"
#endif



double gamma = 1.4;   /* Ratio of specific heats; Adiabatic exponent... */

double* x = 0;

cstate* flux = 0;
cstate* u_old = 0;
cstate* u_new = 0;
pstate* w_new = 0;
pstate* w_old = 0;
pstate* w_intercell = 0; /* Index convention:  W[i] corresponds to W[i-1/2] */
cstate* u_intercell = 0; /* needed for HLL solver */


double t = 0;
double dt = 0;
double dx = 1;

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
  pstate left, right;
  read_ic(&left, &right);



#ifdef RIEMANN_EXACT
  printf("Using exact Riemann solver.\n");
#elif defined RIEMANN_TRRS
  printf("Using TRRS Riemann solver.\n");
#elif defined RIEMANN_TSRS
  printf("Using TSRS Riemann solver.\n");
#elif defined RIEMANN_HLL
  printf("Using HLL Riemann solver.\n");
#elif defined RIEMANN_HLLC
  printf("Using HLLC Riemann solver.\n");
#endif



  /*----------------------------------*/
  /* allocate memory for state arrays */
  /* leave extra space for boundaries */
  /*----------------------------------*/

  x =           malloc((pars.nx+NBCT)*sizeof(double));
  w_old =       malloc((pars.nx+NBCT)*sizeof(pstate));
  w_new =       malloc((pars.nx+NBCT)*sizeof(pstate));
  w_intercell = malloc((pars.nx+NBCT)*sizeof(pstate));
  u_intercell = malloc((pars.nx+NBCT)*sizeof(cstate));
  u_old =       malloc((pars.nx+NBCT)*sizeof(cstate));
  u_new =       malloc((pars.nx+NBCT)*sizeof(cstate));
  flux =        malloc((pars.nx+NBCT)*sizeof(cstate));

  /*----------------------------------*/
  /* initialize values                */
  /*----------------------------------*/

  dx = 2./((double) pars.nx);
  x[0] = -1-NBC*dx;
  for (int i=1; i<pars.nx+NBCT; i++){
    x[i] = x[i-1]+dx;
  }
  double boundary = 0.0; /* separating left and right initial states */

  /* use these in case you want to mimic Toro's E1GODS solver from NUMERICA */
  /* dx = 1./((double) pars.nx); */
  /* x[0] = -NBC*dx; */
  /* for (int i=1; i<pars.nx+NBCT; i++){ */
  /*   x[i] = x[i-1]+dx; */
  /* } */
  /* double boundary = 0.3;  [> separating left and right initial states <] */
 

  if (pars.twopstate_ic){
    for (int i=0; i<pars.nx+NBCT; i++){
      if (x[i]<=boundary) {
        w_old[i].rho = left.rho;
        w_old[i].u   = left.u;
        w_old[i].p   = left.p;
      }
      else {
        w_old[i].rho = right.rho;
        w_old[i].u   = right.u;
        w_old[i].p   = right.p;
      }
    }
  }
  else{
    printf("Can't handle non-Riemann initial conditions yet...\n");
    exit(0);
  }


  /* Initialize counters */
  int step = 0;
  int outputstep = 0;
  int outcount = 0;

  if (pars.verbose) printf("Writing initial output\n");
  write_output(outcount, t, x, w_old);

  compute_conserved_states();

  /* -------------------- */
  /*   Main loop          */
  /* -------------------- */
  while(t < pars.tmax){

    if (pars.nsteps>0 && step == pars.nsteps) break;
    step += 1;
    outputstep += 1;

    compute_fluxes();
    dt=compute_dt(dx);
/* TODO remove */
/* dt = 0.01; */
    if (t+dt > pars.tmax) dt = pars.tmax-t;

    compute_new_states();

/* TODO remove */
/* printf("U_OLD\n"); */
/* for (int i=NBC; i<pars.nx+NBCT; i++){ */
/*   printf("%lf \t %lf \t %lf \n", u_old[i].rho, u_old[i].rhou, u_old[i].E); */
/* } */
/* printf("\n===================================================================\n"); */
/* printf("U_NEW\n"); */
/* for (int i=NBC; i<pars.nx+NBCT; i++){ */
/*   printf("%lf \t %lf \t %lf \n", u_new[i].rho, u_new[i].rhou, u_new[i].E); */
/* } */
/* printf("\n===================================================================\n"); */
/* printf("FLUXES\n"); */
/* for (int i=NBC; i<pars.nx+NBCT-1; i++){ */
/*   printf("%lf \t %lf \t %lf \t %lf \t %lf \t %lf \n", flux[i].rho, flux[i].rhou, flux[i].E, flux[i+1].rho, flux[i+1].rhou, flux[i+1].E); */
/* } */



    /* swap new states with the old ones */
    pstate *ptemp = w_old;
    w_old = w_new;
    w_new = ptemp;

    cstate *ctemp = u_old;
    u_old = u_new;
    u_new = ctemp;

    /* for (int i=NBC; i<pars.nx+NBCT; i++){ */
    /*   w_new[i].rho = 0; */
    /*   w_new[i].u = 0; */
    /*   w_new[i].p = 0; */
    /*   u_new[i].rho = 0; */
    /*   u_new[i].rhou = 0; */
    /*   u_new[i].E = 0; */
    /* } */

    set_boundaries();


    t += dt;
    if (pars.verbose) printf("Finished step %d at t = %10.6lf    dt = %10.6lf\n", step, t, dt);

    if (pars.foutput>0){
      if (outputstep == pars.foutput){
        outputstep = 0;
        outcount += 1;
        if (pars.verbose) printf("Writing output\n");
        write_output(outcount, t, x, w_old);
      }
    }
  }


  if (outputstep>0){
    outcount += 1;
    if (pars.verbose) printf("Writing final output\n");
    write_output(outcount, t, x, w_old);
  }


  printf("Finished clean. Yay!\n");
  return(0);

}




