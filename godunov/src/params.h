/* All around parameters used in the simulation. */

#ifndef PARAMS_H
#define PARAMS_H
#include "gas.h"

#define MAX_FNAME_SIZE 100

/* GLOBAL PARAMETERS */
typedef struct {
  int verbose;                        /* how talkative I am                   */
  int nsteps;                         /* how many steps to take               */
  double tmax;                        /* at what time to end                  */
  int nx;                             /* number of mesh points                */
  double ccfl;                        /* CFL coefficient                      */
  int twostate_ic;                    /* whether IC are left/right state only */
  char paramfilename[MAX_FNAME_SIZE]; /* parameter filename                   */
  char datafilename[MAX_FNAME_SIZE];  /* IC data filename                     */
} params;



// extern void check_params();
extern void init_params();
extern void print_params();


#endif
