#pragma once

#include "vec.h"

vec FindHessEigenValues( matr *m, double step, double eps, int *iter );
vec InverseIter( matr *m, vec *y0, double eps, int *iter );
vec LUMethod( matr *m, double eps, int *iter );
