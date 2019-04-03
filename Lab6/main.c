#include <conio.h>

#include "vec.h"
#include "eigen.h"
#include "decompos.h"

#define SIZE 5
#define STEP 1e-2
#define EPS 1e-8

void main( void )
{
  matr a, a1, rand, rand_inv;
  vec eig, x0 = FullRandVec(SIZE);
  int i, iter;
  double *eigens, eps;

  _CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );

  ///////////// Good separable eigenvalues /////////////
  eigens = malloc(SIZE * sizeof(double));
  for (i = SIZE; i > 0; i--)
    eigens[SIZE - i] = i;

  a = Diag(SIZE, eigens);

  /* random matrix for transformations */
  rand = RandMatr(SIZE);
  rand_inv = Ortogonize(&rand);//InverseMatr(&rand);
  rand = InverseMatr(&rand_inv);

  PrintMatr(&a, "Start diagonal matrix with good separable eigenvalues", stdout);
  //PrintMatr(&rand_inv, "Matrix of eigenvectors", stdout);

  /* R^-1 * A * R */
  a1 = MatrMulMatr(&rand, &a, 0);
  FreeMatr(&a);
  a = MatrMulMatr(&a1, &rand_inv, 0);

  PrintMatr(&a, "Final non-symmetric matrix with same eigenvalues", stdout);

  printf("\nDeterminant: %lf\n", Det(&a));

  PrintVec(&x0, "Start random vector", stdout);

  eig = InverseIter(&a, &x0, EPS, &iter);
  PrintVec(&eig, "First and second minimum eigenvalues", stdout);

  printf("Current epsylon: %g\n", EPS);
  printf("Iterations: %i\n", iter);

  ///////////// Table /////////////
  eps = 0.1;
  printf("\nEpsylon    | Iterations\n");
  printf("-----------------------\n");
  while (eps >= EPS)
  {
    FreeVec(&eig);
    eig = InverseIter(&a, &x0, eps, &iter);
    printf("%.8lf | %i\n", eps, iter);
    eps *= 0.1;
  }

  _getch();

  FreeVec(&eig);
  FreeVec(&x0);
  FreeMatr(&a);
  FreeMatr(&a1);
  FreeMatr(&rand);
  FreeMatr(&rand_inv);
  free(eigens);

  ///////////// Bad separable eigenvalues /////////////
  eigens = malloc(SIZE * sizeof(double));
  for (i = 0; i < SIZE; i++)
    eigens[i] = 1 + pow(0.5, i + 1);

  x0 = FullRandVec(SIZE);

  a = Diag(SIZE, eigens);

  /* random matrix for transformations */
  rand = RandMatr(SIZE);
  rand_inv = InverseMatr(&rand);

  PrintMatr(&a, "Start diagonal matrix with bad separable eigenvalues", stdout);
  //PrintMatr(&rand_inv, "Matrix of eigenvectors", stdout);

  /* R^-1 * A * R */
  a1 = MatrMulMatr(&rand_inv, &a, 0);
  FreeMatr(&a);
  a = MatrMulMatr(&a1, &rand, 0);

  PrintMatr(&a, "Final non-symmetric matrix with same eigenvalues", stdout);
  printf("\nDeterminant: %lf\n", Det(&a));

  PrintVec(&x0, "Start random vector", stdout);

  eig = InverseIter(&a, &x0, EPS, &iter);
  PrintVec(&eig, "First and second minimum eigenvalues", stdout);

  printf("Current epsylon: %g\n", EPS);
  printf("Iterations: %i\n", iter);

  ///////////// Table /////////////
  eps = 0.1;
  printf("\nEpsylon    | Iterations\n");
  printf("-----------------------\n");
  while (eps >= EPS)
  {
    FreeVec(&eig);
    eig = InverseIter(&a, &x0, eps, &iter);
    printf("%.8lf | %i\n", eps, iter);
    eps *= 0.1;
  }

  _getch();

  FreeVec(&eig);
  FreeVec(&x0);
  FreeMatr(&a);
  FreeMatr(&a1);
  FreeMatr(&rand);
  FreeMatr(&rand_inv);
  free(eigens);
}
