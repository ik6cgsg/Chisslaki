#ifndef __VEC_H__
#define __VEC_H__

#pragma once

#include <stdio.h>
#include <crtdbg.h>
#include <stdlib.h>
#include <math.h>

#include "matr.h"

#define RAND_VALUE(MAX) ((rand() * 1.0 / RAND_MAX) * ((MAX) * 2) - (MAX))
#define MAX_RAND_VALUE 5
#define MAX_RAND_OFFSET 0.01

typedef struct
{
  double *v;
  int length;
} vec;

vec Zeros( int length );
vec RandVec( int length );
void FreeVec( vec *v );

vec VecMulNum( vec *v, double n );
vec VecPlusVec( vec *v1, vec *v2 );
vec VecMinusVec( vec *v1, vec *v2 );
matr VecMulVec( vec *v1, vec *v2 );
double Scalar( vec *v1, vec *v2 );
vec MatrMulVec( matr *m, vec *v );
vec OffsetVec( vec *v );

double NormInf( vec *v );
double NormEucl( vec *v );

void PrintVec( vec *v, char *name, FILE *f );

#endif // __VEC_H__
