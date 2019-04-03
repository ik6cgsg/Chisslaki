#ifndef _DIFF_H__
#define _DIFF_H__
#pragma once

#include "grid.h"

typedef double(*funcptr)(double);
typedef double(*funcptr2)(double, double);

grid RungeKutta(funcptr2 diffEq, double y0, grid * gr, uint numInt);
grid AdamsBashforth(funcptr2 diffEq, double y0, grid * gr, uint numInt);
grid AdamsPredCorr(funcptr2 diffEq, double y0, grid * gr, uint numInt);
grid DiffRun(funcptr p, funcptr q, funcptr f, grid * gr, double stCond[6], uint numInt);

#endif // _DIFF_H__
