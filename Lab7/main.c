#include <conio.h>

#include "vec.h"
#include "eigen.h"
#include "decompos.h"

#define SIZE 10
#define STEP 1e-2
#define EPS 1e-8

void main( void )
{
  matr a, a1, rand, rand_inv;
  vec eig;
  int i, iter, j;
  double *eigens, eps;

  _CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );

  ///////////// Good separable eigenvalues /////////////
  eigens = malloc(SIZE * sizeof(double));
  for (i = SIZE; i > 0; i--)
    eigens[SIZE - i] = i;

  a = Diag(SIZE, eigens);

  /* random matrix for transformations */
  rand = RandMatr(SIZE);
  rand_inv = InverseMatr(&rand);

  PrintMatr(&a, "Start diagonal matrix with good separable eigenvalues", stdout);

  /* R^-1 * A * R */
  a1 = MatrMulMatr(&rand, &a, 0);
  FreeMatr(&a);
  a = MatrMulMatr(&a1, &rand_inv, 0);

  PrintMatr(&a, "Final non-symmetric matrix with same eigenvalues", stdout);

  printf("\nDeterminant: %lf\n", Det(&a));

  eig = LUMethod(&a, EPS, &iter);
  PrintVec(&eig, "Eigenvalues", stdout);

  printf("Current epsylon: %g\n", EPS);
  printf("Iterations: %i\n", iter);

  ///////////// Table /////////////
  eps = 0.1;
  printf("\nEpsylon    | Iterations\n");
  printf("-----------------------\n");
  while (eps >= EPS)
  {
    FreeVec(&eig);
    eig = LUMethod(&a, eps, &iter);
    printf("%.8lf | %i\n", eps, iter);
    eps *= 0.1;
  }

  _getch();

  FreeVec(&eig);
  FreeMatr(&a);
  FreeMatr(&a1);
  FreeMatr(&rand);
  FreeMatr(&rand_inv);
  free(eigens);

  ///////////// Bad separable eigenvalues /////////////
  eigens = malloc(SIZE * sizeof(double));
  /*
  for (i = 0; i < SIZE; i++)
  {
    //eigens[i] = 1 + pow(0.3, i + 1);

    if (i == SIZE / 2)
      eigens[i] = i + 0.00001;//1 + pow(0.3, i + 1);
    else
      eigens[i] = i + 1;
    
  } */
  for (i = SIZE; i > 0; i--)
  {
    if (i == SIZE / 2)
      eigens[SIZE - i] = (i - 1) + 0.003;
    else
      eigens[SIZE - i] = i;
  }

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

  eig = LUMethod(&a, EPS, &iter);
  PrintVec(&eig, "Eigenvalues", stdout);

  printf("Current epsylon: %g\n", EPS);
  printf("Iterations: %i\n", iter);

  ///////////// Table /////////////
  eps = 0.1;
  printf("\nEpsylon    | Iterations\n");
  printf("-----------------------\n");
  while (eps >= EPS)
  {
    FreeVec(&eig);
    eig = LUMethod(&a, eps, &iter);
    printf("%.8lf | %i\n", eps, iter);
    eps *= 0.1;
  }

  _getch();

  FreeVec(&eig);
  FreeMatr(&a);
  FreeMatr(&a1);
  FreeMatr(&rand);
  FreeMatr(&rand_inv);
  free(eigens);
}
