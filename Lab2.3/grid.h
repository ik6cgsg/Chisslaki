#ifndef _GRID_H__
#define _GRID_H__
#pragma once

#include <math.h>

typedef unsigned int uint;

typedef enum
{
  GRID_UNIFORM, GRID_CHEBYSHEV, GRID_USER
} grid_type; /* end of 'grid_type' enum */

typedef struct
{
  uint size;
  double *nodes;
} grid; /* end of 'grid' struct */

int GridInit( grid *gr, uint size );

void GridUniform( grid *gr, double begin, double end );
void GridChebyshev( grid *gr, double begin, double end );
void GridHilb( grid *gr, double begin, double end );

void GridFree( grid *gr );

#endif /* _GRID_H__ */
