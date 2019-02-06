#ifndef GODUNOV_H
#define GODUNOV_H

extern void compute_intercell_states();
extern void compute_fluxes();
extern void compute_new_states();
extern void set_boundaries();
extern void compute_conserved_states();

#endif
