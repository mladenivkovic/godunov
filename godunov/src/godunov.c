#include <stdio.h>

#include "gas.h"
#include "godunov.h"
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

extern pstate* w_old;
extern pstate* w_new;
extern pstate* w_intercell;
extern cstate* u_old;
extern cstate* u_new;
extern cstate* u_intercell;
extern cstate* flux;
extern params pars;
extern double dx;
extern double dt;
extern double gamma;

extern double *x;





/* ===================================== */
void compute_intercell_states(){
/* ===================================== */
  /* Compute intercell states, based on  */
  /* whether we have vacuum or not       */
  /* first compute primitives, then      */
  /* convert to conserved variables      */
  /*-------------------------------------*/

  pstate left, right, starL, starR;

  for (int i = 2; i<pars.nx+3; i++){
    left = w_old[i-1];
    right = w_old[i];
    int vacuum = check_vacuum(&left, &right);

    if (vacuum){
      compute_riemann_vacuum(&left, &right, &w_intercell[i]);
    }
    else {
      compute_star_pstate(&left, &right, &starL, &starR);
      compute_riemann(&left, &right, &starL, &starR, &w_intercell[i]);
    }

    u_intercell[i].rho = w_intercell[i].rho;
    u_intercell[i].rhou = w_intercell[i].rho * w_intercell[i].u;
    u_intercell[i].E = energy(&w_intercell[i]);
  }

}




/* =================================== */
void compute_fluxes(){
/* =================================== */
  /* compute the intercell fluxes      */
  /* fluxes are for conserved          */
  /* variables, not primitives!!!!     */
  /*-----------------------------------*/

  for (int i=1; i<pars.nx+3; i++){
    pstate ic = w_intercell[i];

    flux[i].rho = ic.rho * ic.u;
    flux[i].rhou = ic.rho * ic.u*ic.u + ic.p;
    flux[i].E = ic.u*(energy(&ic) + ic.p);

  }
}



/* ================================ */
void compute_conserved_states(){
/* ================================ */
  /* Computes conserved states from */
  /* primitive states               */
  /*--------------------------------*/

  for (int i = 0; i<pars.nx+4; i++){
    u_old[i].rho = w_old[i].rho;
    u_old[i].rhou = w_old[i].rho * w_old[i].u;
    u_old[i].E = energy(&w_old[i]);
  }


}

/* =================================== */
void compute_new_states(){
/* =================================== */
  /* Compute the new pstates           */
  /*-----------------------------------*/

  double dtdx = dt/dx;

  for (int i=2; i<pars.nx+2; i++){
    u_new[i].rho = u_old[i].rho + dtdx*(flux[i].rho - flux[i+1].rho);
    u_new[i].rhou = u_old[i].rhou + dtdx*(flux[i].rhou - flux[i+1].rhou);
    u_new[i].E = u_old[i].E + dtdx*(flux[i].E - flux[i+1].E);

    w_new[i].rho = u_new[i].rho;
    w_new[i].u = u_new[i].rhou/u_new[i].rho;
    w_new[i].p = (gamma-1)*(u_new[i].E - 0.5*u_new[i].rhou*w_new[i].u);
  }

}


/* ==========================*/
void set_boundaries(){
/* ==========================*/
  /* Set boundary conditions */
  /*-------------------------*/

  /* transmissive boundaries */

  w_old[0].rho = w_old[3].rho;
  w_old[0].u = w_old[3].u;
  w_old[0].p = w_old[3].p;

  w_old[1].rho = w_old[2].rho;
  w_old[1].u = w_old[2].u;
  w_old[1].p = w_old[2].p;

  w_old[pars.nx+2].rho = w_old[pars.nx+1].rho;
  w_old[pars.nx+2].u = w_old[pars.nx+1].u;
  w_old[pars.nx+2].p = w_old[pars.nx+1].p;

  w_old[pars.nx+3].rho = w_old[pars.nx].rho;
  w_old[pars.nx+3].u = w_old[pars.nx].u;
  w_old[pars.nx+3].p = w_old[pars.nx].p;

}
