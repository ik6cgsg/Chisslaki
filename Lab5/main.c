#include <conio.h>

#include "vec.h"
#include "eigen.h"
#include "decompos.h"

#define SIZE 5
#define STEP 1e-2
#define EPS 1e-7

void main( void )
{
  matr a, a1, rand, rand_inv, h;
  vec eig;
  int i, iter;
  double *eigens;

  _CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );

  ///////////// Good separable eigenvalues /////////////
  eigens = malloc(SIZE * sizeof(double));
  for (i = 0; i < SIZE; i++)
    eigens[i] = i + 1;

  a = Diag(SIZE, eigens);

  /* random matrix for transformations */
  rand = RandMatr(SIZE);
  rand_inv = InverseMatr(&rand);

  PrintMatr(&a, "Start diagonal matrix with good separable eigenvalues", stdout);

  /* R^-1 * A * R */
  a1 = MatrMulMatr(&rand_inv, &a, 0);
  FreeMatr(&a);
  a = MatrMulMatr(&a1, &rand, 0);

  PrintMatr(&a, "Final non-symmetric matrix with same eigenvalues", stdout);
  printf("\nDeterminant: %lf\n", Det(&a));

  h = ToHess(&a);
  PrintMatr(&h, "Hessenberg", stdout);

  eig = FindHessEigenValues(&a, STEP, EPS, &iter);
  PrintVec(&eig, "Eigenvalues of Hessenberg", stdout);
  printf("\nIterations: %i\n", iter);

  _getch();

  FreeVec(&eig);
  FreeMatr(&a);
  FreeMatr(&a1);
  FreeMatr(&rand);
  FreeMatr(&rand_inv);
  FreeMatr(&h);
  free(eigens);

  ///////////// Bad separable eigenvalues /////////////
  eigens = malloc(SIZE * sizeof(double));
  for (i = 0; i < SIZE; i++)
    eigens[i] = 1 + pow(1e-1, i + 1);

  a = Diag(SIZE, eigens);

  /* random matrix for transformations */
  rand = RandMatr(SIZE);
  rand_inv = InverseMatr(&rand);

  PrintMatr(&a, "Start diagonal matrix with bad separable eigenvalues", stdout);

  /* R^-1 * A * R */
  a1 = MatrMulMatr(&rand_inv, &a, 0);
  FreeMatr(&a);
  a = MatrMulMatr(&a1, &rand, 0);

  PrintMatr(&a, "Final non-symmetric matrix with same eigenvalues", stdout);
  printf("\nDeterminant: %lf\n", Det(&a));

  h = ToHess(&a);
  PrintMatr(&h, "Hessenberg", stdout);

  eig = FindHessEigenValues(&a, 1e-4, EPS, &iter);
  PrintVec(&eig, "Eigenvalues of Hessenberg", stdout);
  printf("\nIterations: %i\n", iter);

  _getch();

  FreeVec(&eig);
  FreeMatr(&a);
  FreeMatr(&a1);
  FreeMatr(&rand);
  FreeMatr(&rand_inv);
  FreeMatr(&h);
  free(eigens);
}
