#ifndef __MATR_H__
#define __MATR_H__

#pragma once

#include <stdio.h>

typedef struct
{
  double **a;
  int size;
} matr;

matr Unit( int size );
matr RandMatr( int size );
matr Hilb( int size );
void FreeMatr( matr *m );

matr MakeMatr( int size, double *a );
matr MatrMulNum( matr *m, double n );
matr MatrMulMatr( matr *m1, matr *m2, int ind );
matr MatrPlusMatr( matr *m1, matr *m2 );
matr MatrMinusMatr( matr *m1, matr *m2 );
matr InverseMatr( matr *m );
matr OffsetMatr( matr *m );
matr TransposeMatr( matr *m );
matr MakeSymPosMatr( int size );
void CopyMatr( matr *m1, matr *m2 );

void PrintMatr( matr *m, char *name, FILE *f );

void SwapPtr( double **a, double **b );

double Det( matr *m );
double GetCond( matr *m );
double InfNorm( matr *m );

#endif // __MATR_H__
