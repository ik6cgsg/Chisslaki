#ifndef _INTERPOL_H__
#define _INTERPOL_H__
#pragma once

#include "grid.h"

double InterPolLagrange( double (*f)(double), grid const *gr, double x );
void GenerateY( grid const *x, grid *y, double (*f)(double), uint nodes_num, grid_type gt);

#endif /* _INTERPOL_H__ */
