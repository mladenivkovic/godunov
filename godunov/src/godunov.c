#include "gas.h"
#include "godunov.h"
#include "params.h"
#if RIEMANN==EXACT
#include "riemann-exact.h"
#endif

extern state* w_old;
extern state* w_new;
extern state* w_intercell;
extern params pars;
extern double dx;
extern double vmax;






void compute_intercell_states(){

  state left, right, starL, starR;

  for (int i = 1; i<pars.nx+3; i++){
    left = w_old[i];
    right = w_old[i+1];
    int vacuum = check_vacuum(&left, &right);

    if (vacuum){
      w_intercell[i] = compute_riemann_vacuum(&left, &right);
    }
    else {
      compute_star_state(&left, &right, &starL, &starR);
      w_intercell[i] = compute_riemann(&left, &right, &starL, &starR);
    }

  }

}
