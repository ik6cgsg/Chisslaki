#ifndef _DIFF_H__
#define _DIFF_H__
#pragma once

#include "grid.h"

grid RungeKutta(double(*diffEq)(double, double), double y0, grid * gr, uint numInt);
grid AdamsBashforth(double(*diffEq)(double, double), double y0, grid * gr, uint numInt);
grid AdamsPredCorr(double(*diffEq)(double, double), double y0, grid * gr, uint numInt);

#endif // _DIFF_H__
