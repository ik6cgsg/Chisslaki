#include <stdio.h>
#include <math.h>

#include "integr.h"
#include "func.h"

void main( void )
{
  int i;
  double Int1, Int2;

  printf("first function: e^x for x=(-1, 1)\n\n");

  printf("intervals | 3/8 result | Newton-Leibniz | practic error | theory error\n");
  printf("----------------------------------------------------------------------\n");

  Int1 = NewLeyb(func1int, -1, 1);

  for (i = 4; i <= 80; i++)
  {
    double res = Integrate(func1, -1, 1, i);
    double h = 2.0 / (i - 1) / 3;

    if (i > 20)
      i += 5;

    printf("%3i       | %2.5lf    | %g         | %e  | %g\n", i - 1, res, Int1, fabs(Int1 - res), 3 * pow(h, 5) / 80 * func1deriv4(1));
  }

  printf("\nsecond function: (e^x)^5 for x=(-1, 1)\n\n");

  printf("intervals | 3/8 result | Newton-Leibniz | practic error | theory error\n");
  printf("----------------------------------------------------------------------\n");

  Int2 = NewLeyb(func2int, -1, 1);

  for (i = 4; i <= 80; i++)
  {
    double res = Integrate(func2, -1, 1, i);
    double h = 2.0 / (i - 1) / 3;

    if (i > 20)
      i += 5;

    printf("%3i       | %2.5lf   | %g        | %e  | %g\n", i - 1, res, Int2, fabs(Int2 - res), 3 * pow(h, 5) / 80 * func2deriv4(1));
  }

  printf("\n\nCounting by given epsylon\n\n");

  printf("first function: e^x for x=(-1, 1)\n\n");

  printf("  epsylon   | 3/8 result | intervals\n");
  printf("-----------------------------------\n");

  for (i = 1; i <= 10; i++)
  {
    double eps = pow(0.1, i);
    int N;
    double res = IntegrateRunge(func1, -1, 1, eps, &N);

    printf("%e| %lf   | %i\n", eps, res, N - 1);
  }

  printf("\nsecond function: (e^x)^5 for x=(-1, 1)\n\n");

  printf("  epsylon   | 3/8 result | intervals\n");
  printf("-----------------------------------\n");

  for (i = 1; i <= 10; i++)
  {
    double eps = pow(0.1, i);
    int N;
    double res = IntegrateRunge(func2, -1, 1, eps, &N);

    printf("%e| %lf  | %i\n", eps, res, N - 1);
  }

  getchar();
}
