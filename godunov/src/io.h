/* IO routines */

#ifndef IO_H
#define IO_H
#include "gas.h"
#include "params.h"

extern void read_cmdlineargs(int argc, char* argv[]);
extern void read_ic(state* left, state* right);
extern void read_paramfile();
extern void write_output(int step, double t, double* x, double* rho, double* u, double* p);

#endif
