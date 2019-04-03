#include <stdio.h>

#include "rootf.h"

#define TRESHOLD 1e-10
#define ITER_NUM 6

double PolyFun( double x )
{
  return x * x * x - 2 * x * x - 4 * x + 7;
}

double PolyPhi1( double x )
{
  return x - 1 / 14.7 * PolyFun(x);
}

double PolyPhi2( double x )
{
  return x + 1 / 3.6 * PolyFun(x);
}

double PolyPhi3( double x )
{
  return x - 1 / 4.7 * PolyFun(x);
}

double TransFun( double x )
{
  return 2 * x * sin(x) - cos(x);
}

double TransPhi1( double x )
{
  return x + 2 / 5.7 * TransFun(x);
}

double TransPhi2( double x )
{
  return x - 2 / 5.7 * TransFun(x);
}

double TransPhi3( double x )
{
  return x + 2 / 14.0 * TransFun(x);
}

int main( void )
{
  double x;
  int n = 0, i;
  double eps = 0.1;

  /* Plynomial */
  printf("POLYNOM\n\n\n");

  /* Half divide mehtod */
  printf("half division\n\n");
  if (n = HalfDiv(PolyFun, -3, 0, &x, TRESHOLD))
    printf("root1 is %.10lf  iterations: %i  epsylon: %g\n", x, n, TRESHOLD);
  if (n = HalfDiv(PolyFun, -2 / 3.0, 2, &x, TRESHOLD))
    printf("root2 is %.10lf  iterations: %i  epsylon: %g\n", x, n, TRESHOLD);
  if (n = HalfDiv(PolyFun, 2.1, 10, &x, TRESHOLD))
    printf("root3 is %.10lf  iterations: %i  epsylon: %g\n", x, n, TRESHOLD);

  /* Simple iterations method */
  printf("simple iterations\n\n");
  if (n = SimpleIter(PolyPhi1, -2, 0.029, &x, TRESHOLD))
    printf("root1 is %.10lf  iterations: %i  epsylon: %g\n", x, n, TRESHOLD);
  if (n = SimpleIter(PolyPhi2, 1.5, 0.059, &x, TRESHOLD))
    printf("root2 is %.10lf  iterations: %i  epsylon: %g\n", x, n, TRESHOLD);
  if (n = SimpleIter(PolyPhi3, 2.5, 0.075, &x, TRESHOLD))
    printf("root3 is %.10lf  iterations: %i  epsylon: %g\n", x, n, TRESHOLD);

   /* Transcendent */
  printf("TRANSCENDENT\n\n\n");

  /* Half divide mehtod */
  printf("half division\n\n");
  if (n = HalfDiv(TransFun, -1.5, 0, &x, TRESHOLD))
    printf("root1 is %.10lf  iterations: %i  epsylon: %g\n", x, n, TRESHOLD);
  if (n = HalfDiv(TransFun, 0.4, 1, &x, TRESHOLD))
    printf("root2 is %.10lf  iterations: %i  epsylon: %g\n", x, n, TRESHOLD);
  if (n = HalfDiv(TransFun, 2, 4, &x, TRESHOLD))
    printf("root3 is %.10lf  iterations: %i  epsylon: %g\n", x, n, TRESHOLD);

  /* Simple iterations method */
  printf("simple iterations\n\n");
  if (n = SimpleIter(TransPhi1, -2, 0.015, &x, TRESHOLD))
    printf("root1 is %.10lf  iterations: %i  epsylon: %g\n", x, n, TRESHOLD);
  if (n = SimpleIter(TransPhi2, 1, 0.015, &x, TRESHOLD))
    printf("root2 is %.10lf  iterations: %i  epsylon: %g\n", x, n, TRESHOLD);
  if (n = SimpleIter(TransPhi3, 3.3, 0.011, &x, TRESHOLD))
    printf("root3 is %.10lf  iterations: %i  epsylon: %g\n\n\n\n", x, n, TRESHOLD);

  /* Treshold table */
  printf("\nPOLYNOM\nHalf division\n");
  for (i = 1; i < ITER_NUM; i++)
  {
    eps *= 0.01;
    printf("root is %.10lf  iterations: %i  epsylon: %g\n", x, HalfDiv(PolyFun, -3, 0, &x, eps), eps);
  }
  eps = 0.1;
  printf("Simple iterations\n");
  for (i = 1; i < ITER_NUM; i++)
  {
    eps *= 0.01;
    printf("root is %.10lf  iterations: %i  epsylon: %g\n", x, SimpleIter(PolyPhi1, -2, 0.029, &x, eps), eps);
  }
  eps = 0.1;
  printf("\nTRANSCENDENT\nSimple iterations\n");
  for (i = 1; i < ITER_NUM; i++)
  {
    eps *= 0.01;
    printf("root is %.10lf  iterations: %i  epsylon: %g\n", x, HalfDiv(TransFun, 0.4, 1, &x, eps), eps);
  }
  eps = 0.1;
  printf("Simple iterations\n");
  for (i = 1; i < ITER_NUM; i++)
  {
    eps *= 0.01;
    printf("root is %.10lf  iterations: %i  epsylon: %g\n", x, SimpleIter(TransPhi2, 1, 0.015, &x, eps), eps);
  }

  return 0;
}
