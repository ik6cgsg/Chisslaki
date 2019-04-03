#ifndef _FUNC_H__
#define _FUNC_H__
#pragma once

#include "spline.h"

double func1( double x );
double func2( double x );
double secDerivFunc2( double x );
double maxDelta( double (*f)(double), spline *sp );
void generateFile( char const *name, int points, int nodes );

#endif /* _FUNC_H__ */
