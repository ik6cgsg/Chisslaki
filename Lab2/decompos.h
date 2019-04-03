#ifndef __DECOMPOS_H__
#define __DECOMPOS_H__

#pragma once

#include "matr.h"
#include "vec.h"

int LUDecompose( matr *orig, matr *l, matr *u );
int QRDecompose( matr *orig, matr *q, matr*r );
vec SolveLU( matr *m, vec *b );
vec SolveQR( matr *m, vec *b );

#endif;
