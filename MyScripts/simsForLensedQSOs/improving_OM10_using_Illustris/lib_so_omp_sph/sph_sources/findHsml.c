#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "proto.h"
#include "allvars_SPH.h"
#include "intfuncs.h"

int findHsml(PARTICLE *particle, long *NumP, long * NgbP, double *BS, float * SmoothLength)
{
  NumPart = * NumP;
  P = (struct particle_data *)particle;
  DesDensNgb = * NgbP;
  BoxSize = * BS;
  Hsml = SmoothLength;
  peano_hilbert_order();
  tree_treeallocate(2.0 * NumPart, NumPart);
  Softening = 0.05;
  tree_treebuild();
  determine_hsml();
  tree_treefree();
  return 0;
}


#ifdef PERIODIC
#define NGB_PERIODIC(x) (((x)>BoxHalf)?((x)-BoxSize):(((x)<-BoxHalf)?((x)+BoxSize):(x)))
#else
#define NGB_PERIODIC(x) (x)
#endif
