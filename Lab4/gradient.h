#pragma once

#include "vec.h"

vec Grad( matr * a, vec * b, vec * x );
vec SolveGrad( matr * a, vec * b, vec * x0, double eps, int *iter );
