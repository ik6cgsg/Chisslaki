#include <conio.h>

#include "vec.h"
#include "gradient.h"

#define SIZE 10
#define EPS 1e-6

void main( void )
{
  matr a;
  vec b, x, x0 = RandVec(SIZE), x_sol = RandVec(SIZE), ax;
  int iter, i;
  double eps;

  _CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );

  PrintVec(&x_sol, "Given solution", stdout);

  printf("___________________________________________________________");

  a = MakeSymPosMatr(SIZE);
  b = MatrMulVec(&a, &x_sol);

  PrintMatr(&a, "Random symmetric positive matrix", stdout);
  printf("\nDet is %g", Det(&a));
  printf("\nCond is %g", GetCond(&a));
  printf("\nEpsylon is %g", EPS);

  PrintVec(&b, "\nRight part", stdout);

  x = SolveGrad(&a, &b, &x0, EPS, &iter);
  PrintVec(&x, "solution", stdout);
  printf("iterations: %i\n", iter);

  ax = MatrMulVec(&a, &x);
  FreeVec(&x);
  x = VecMinusVec(&ax, &b);
  FreeVec(&ax);

  PrintVec(&x, "Discrepancy", stdout);

  //_getch();

  FreeMatr(&a);
  FreeVec(&b);
  FreeVec(&x);
  //FreeVec(&x0);
  //FreeVec(&x_sol);

  printf("___________________________________________________________");

  a = Hilb(SIZE);
  b = MatrMulVec(&a, &x_sol);

  PrintMatr(&a, "Hilbert matrix", stdout);
  printf("\nDet is %g", Det(&a));
  printf("\nCond is %g", GetCond(&a));
  printf("\nEpsylon is %g", EPS);

  PrintVec(&b, "\nRight part", stdout);

  x = SolveGrad(&a, &b, &x0, EPS, &iter);
  PrintVec(&x, "solution", stdout);
  printf("iterations: %i\n", iter);

  ax = MatrMulVec(&a, &x);
  FreeVec(&x);
  x = VecMinusVec(&ax, &b);
  FreeVec(&ax);

  PrintVec(&x, "Discrepancy", stdout);

  _getch();

  FreeMatr(&a);
  FreeVec(&b);
  FreeVec(&x);
  //FreeVec(&x0);
  //FreeVec(&x_sol);

  printf("___________________________________________________________");

  a = MakeSymPosMatr(SIZE);
  b = MatrMulVec(&a, &x_sol);

  eps = 1e-2;

  printf("\nIterations table for first matrix:\nEpsylon  | Iterations\n---------------------\n");
  for (i = 0; i < 10; i++)
  {
    eps /= 2;
    x = SolveGrad(&a, &b, &x0, eps, &iter);
    FreeVec(&x);

    printf("%lf | %i\n", eps, iter);
  }
  printf("---------------------\n");

  FreeMatr(&a);
  FreeVec(&b);
  //FreeVec(&x);
  //FreeVec(&x0);
  //FreeVec(&x_sol);

  a = Hilb(SIZE);
  b = MatrMulVec(&a, &x_sol);

  eps = 1e-2;

  printf("\nIterations table for hilbert matrix:\nEpsylon  | Iterations\n---------------------\n");
  for (i = 0; i < 10; i++)
  {
    eps /= 2;
    x = SolveGrad(&a, &b, &x0, eps, &iter);
    FreeVec(&x);

    printf("%lf | %i\n", eps, iter);
  }
  printf("---------------------\n");

  _getch();

  FreeMatr(&a);
  FreeVec(&b);
  FreeVec(&x0);
  FreeVec(&x_sol);
}
