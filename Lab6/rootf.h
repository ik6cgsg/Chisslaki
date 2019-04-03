#ifndef _ROOTF_H_INCLUDED__
#define _ROOTF_H_INCLUDED__

#pragma once

#include <math.h>

int HalfDiv( double (*f)(double), double left, double right, double *x, double TRESHOLD );
int SimpleIter( double (*f)(double), double x0, double q, double *x, double TRESHOLD );

#endif // _ROOTF_H_INCLUDED__
