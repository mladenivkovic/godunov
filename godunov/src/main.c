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
pstate* w_new = 0;
pstate* w_old = 0;
pstate* w_intercell = 0; /* W[i] corresponds to W[i-1/2] */
cstate* u_old = 0;
cstate* u_new = 0;
cstate* u_intercell = 0;
cstate* flux = 0;

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
  pstate left, right;
  read_ic(&left, &right);



  /* allocate memory for pstate arrays */
  /* leave extra space for boundaries */
  /* initialize values                */

  x = malloc((pars.nx+4)*sizeof(double));
  dx = 2./((double) pars.nx);
  x[0] = -1-2*dx;
  for (int i=1; i<pars.nx+4; i++){
    x[i] = x[i-1]+dx;
  }

  w_old = malloc((pars.nx+4)*sizeof(pstate));
  w_new = malloc((pars.nx+4)*sizeof(pstate));
  w_intercell = malloc((pars.nx+4)*sizeof(pstate));
  u_old = malloc((pars.nx+4)*sizeof(cstate));
  u_new = malloc((pars.nx+4)*sizeof(cstate));
  u_intercell = malloc((pars.nx+4)*sizeof(cstate));
  flux = malloc((pars.nx+4)*sizeof(cstate));

  for (int i=0; i<pars.nx+4; i++){
    if (x[i]<0.0) {
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



  int step = 0;
  if (pars.foutput==0) pars.foutput=1;
  int outputstep = 0;
  int outcount = 0;

  if (pars.verbose) printf("Writing initial output\n");
  write_output(outcount, t, x, w_old);

  while(t < pars.tmax){

    if (pars.nsteps>0 && step == pars.nsteps) break;
    step += 1;
    outputstep += 1;

    /* reset vmax */
    vmax = 0;

    compute_conserved_states();
    compute_intercell_states();
    compute_fluxes();

    /* compute new dt */
    dt = pars.ccfl * dx/vmax;

    compute_new_states();

    /* swap new states with the old ones */
    pstate * ptemp = w_old;
    w_old = w_new;
    w_new = ptemp;

    cstate * ctemp = u_old;
    u_old = u_new;
    u_new = ctemp;

    set_boundaries();

    t += dt;
    if (pars.verbose) printf("Finished step %d at t = %10.6lf    dt = %10.6lf\n", step, t, dt);

    if (outputstep == pars.foutput || t >= pars.tmax){
      outputstep = 0;
      outcount += 1;
      if (pars.verbose) printf("Writing output\n");
      write_output(outcount, t, x, w_old);
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




