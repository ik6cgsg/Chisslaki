#ifndef _SPLINE_H__
#define _SPLINE_H__
#pragma once

#include "grid.h"

typedef struct
{
  double a, b, c, d;
} cub_poly; /* end of 'cub_poly' struct */

typedef struct
{
  cub_poly *polys;
  grid gr;
} spline; /* end of 'spline' struct */

void SplineInit( spline *sp, uint size, double begin, double end );

void BuildNaturalSpline( spline *sp, double (*f)(double) );
void BuildSecDerivSpline( spline *sp, double (*f)(double), double sec_deriv_x0, double sec_deriv_xn );
double GetValue( spline *sp, double x );

void SplineFree( spline *sp );

#endif /* _SPLINE_H__ */
