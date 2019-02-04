/* all stuff concerning Riemann problem */

#ifndef RIEMANN_EXACT_H
#define RIEMANN_EXACT_H
#include "gas.h"


extern int check_vacuum(state *left, state *right);
extern void compute_star_state(state *left, state *right, state* starL, state* starR);
extern double compute_pstar(state *left, state *right);
extern double fp(double pguess, state *s, double gamma, double A, double B, double a);
extern double dfpdp(double pguess, state *s, double gamma, double A, double B, double a);
extern double fpfull(double fpL, double fpR, double delta_u);
extern double dfpdpfull(double dfpdpL, double dfpdpR);
extern double rho_star(state *s, state *star);
extern state compute_riemann(state* left, state* right, state* starL, state* starR);
extern state compute_riemann_vacuum(state* left, state* right);
extern double rho_fanL(state* s);
extern double u_fanL(state* s);
extern double p_fanL(state* s);
extern double rho_fanR(state* s);
extern double u_fanR(state* s);
extern double p_fanR(state* s);
#endif
