#include <stdio.h>
#include <math.h>

#include "integr.h"
#include "func.h"

void main( void )
{
  int i;
  double I1 = NewLeib(func1int, -1, 1), I2 = NewLeib(func2int, -1, 1);

  printf("first function: e^x for x=(-1, 1)\n\n");

  printf("intervals |  3/8 result   | Gauss result |   3/8 error  | Gauss error\n");
  printf("-----------------------------------------------------------------------\n");

  for (i = 2; i <= 100; i++)
  {
    double res = Integrate38(func1, -1, 1, i);
    double res2 = IntegrateGauss(func1, -1, 1, i);

    printf("%3i       | %2.8lf    | %2.8lf   | %e | %e\n", i - 1, res, res2, fabs(I1 - res), fabs(I1 - res2));

    if (i > 20)
      i += 5;
  }

  printf("\nsecond function: (e^x)^5 for x=(-1, 1)\n\n");

  printf("intervals |  3/8 result   | Gauss result |   3/8 error  | Gauss error\n");
  printf("-----------------------------------------------------------------------\n");

  for (i = 4; i <= 100; i++)
  {
    double res = Integrate38(func2, -1, 1, i);
    double res2 = IntegrateGauss(func2, -1, 1, i);

    printf("%3i       | %2.8lf   | %2.8lf  | %e | %e \n", i - 1, res, res2, fabs(I2 - res), fabs(I2 - res2));

    if (i > 20)
      i += 5;
  }

  printf("\n\nCounting by given epsylon\n\n");

  printf("first function: e^x for x=(-1, 1)\n\n");

  printf("  epsylon    | 3/8 result | 3/8 intervals | Gauss result | Gauss intervals\n");
  printf("---------------------------------------------------------------------------\n");

  for (i = 1; i <= 10; i++)
  {
    double eps = pow(0.1, i);
    int N, N2;
    double res = Integrate38ByEps(func1, -1, 1, eps, &N);
    double res2 = IntegrateGaussByEps(func1, -1, 1, eps, I1, &N2);

    printf("%e | %2.8lf |     %3i       | %2.8lf   |    %3i\n", eps, res, N - 1, res2, N2 - 1);
  }

  printf("\nsecond function: (e^x)^5 for x=(-1, 1)\n\n");

  printf("  epsylon    | 3/8 result  | 3/8 intervals | Gauss result | Gauss intervals\n");
  printf("---------------------------------------------------------------------------\n");

  for (i = 1; i <= 10; i++)
  {
    double eps = pow(0.1, i);
    int N, N2;
    double res = Integrate38ByEps(func2, -1, 1, eps, &N);
    double res2 = IntegrateGaussByEps(func2, -1, 1, eps, I2, &N2);

    printf("%e | %2.8lf |      %3i      | %2.8lf  |    %3i\n", eps, res, N - 1, res2, N2 - 1);
  }

  getchar();
}
