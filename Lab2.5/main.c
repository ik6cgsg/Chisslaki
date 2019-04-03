#include <stdio.h>
#include <math.h>

#include "diff.h"

double func(double x)
{
  return exp(-3 / x) * x * x;
}

double diffEq(double x, double y)
{
  return (3 + 2 * x) / (x * x) * y;
}

void main(void)
{
  grid grY, grY2, grX;
  uint i, j;

  GridInit(&grX, 5);
  GridUniform(&grX, 1, 4);

  /* Runge-Knutta */
#if 1
  printf("solution by Runge-Kutta (p = 4) for d.e.: x2y' - 2xy = 3y, x = [1, 4]\nintervals - num of intervals between nodes\n\n");
  printf("intervals | node | Runge-Knutta |    exact    | error\n");
  printf("-------------------------------------------------------------\n");
  for (i = 1; i <= 50; i += 7)
  {
    double a, b;
    grY = RungeKutta(diffEq, func(grX.nodes[0]), &grX, i);


    for (j = 0; j < grX.size; j++)
    {
      if (j == grX.size / 2)
        printf("%6i    | %4g | %2.9lf  | %2.9lf | %g\n", i, grX.nodes[j], a = grY.nodes[j], b = func(grX.nodes[j]), fabs(a - b));

      printf("          | %4g | %2.9lf  | %2.9lf | %g\n", grX.nodes[j], a = grY.nodes[j], b = func(grX.nodes[j]), fabs(a - b));
    }

    printf("-------------------------------------------------------------\n");
    GridFree(&grY);
  }


  printf("\n\nerror in start condition = 1%%\n\n");
  printf("intervals | node | Runge-Knutta | without error | difference (%%)\n");
  printf("-------------------------------------------------------------\n");
  for (i = 1; i <= 50; i += 7)
  {
    double a, b;
    grY = RungeKutta(diffEq, func(grX.nodes[0]), &grX, i);
    grY2 = RungeKutta(diffEq, (a = func(grX.nodes[0])) + 0.01 * a, &grX, i);


    for (j = 0; j < grX.size; j++)
    {
      if (j == grX.size / 2)
        printf("%6i    | %4g | %2.9lf  |  %2.9lf  | %g%%\n", i, grX.nodes[j], a = grY2.nodes[j], b = grY.nodes[j], fabs(a - b) / b * 100);

      printf("          | %4g | %2.9lf  |  %2.9lf  | %g%%\n", grX.nodes[j], a = grY2.nodes[j], b = grY.nodes[j], fabs(a - b) / b * 100);
    }

    printf("-------------------------------------------------------------\n");
    GridFree(&grY);
    GridFree(&grY2);
  }

  GridFree(&grX);
#endif
  /* Adams-Bashforth */
#if 0
  printf("solution by Adams (k = 3) for d.e.: x2y' - 2xy = 3y, x = [1, 4]\nintervals - num of intervals between nodes\n\n");
  printf("intervals | node |    Adams     |    exact    | error\n");
  printf("-------------------------------------------------------------\n");
  for (i = 1; i <= 50; i += 7)
  {
    double a, b;
    grY = AdamsBashforth(diffEq, func(grX.nodes[0]), &grX, i);

    for (j = 0; j < grX.size; j++)
      if (i == 1)
        if (j == grX.size / 2)
          printf("     1    | %4g |      --      | %2.9lf |      --      \n", grX.nodes[j], func(grX.nodes[j]));
        else
          printf("          | %4g |      --      | %2.9lf |      --      \n", grX.nodes[j], func(grX.nodes[j]));
      else
      {
        if (j == grX.size / 2)
          printf("%6i    | %4g | %2.9lf  | %2.9lf | %g\n", i, grX.nodes[j], a = grY.nodes[j], b = func(grX.nodes[j]), fabs(a - b));

        printf("          | %4g | %2.9lf  | %2.9lf | %g\n", grX.nodes[j], a = grY.nodes[j], b = func(grX.nodes[j]), fabs(a - b));
      }

    printf("-------------------------------------------------------------\n");
    GridFree(&grY);
  }

  printf("\n\nerror in start condition = 1%%\n\n");
  printf("intervals | node |    Adams     | without error | difference (%%)\n");
  printf("-------------------------------------------------------------\n");
  for (i = 1; i <= 50; i += 7)
  {
    double a, b;
    grY = AdamsBashforth(diffEq, func(grX.nodes[0]), &grX, i);
    grY2 = AdamsBashforth(diffEq, (a = func(grX.nodes[0])) + 0.01 * a, &grX, i);

    for (j = 0; j < grX.size; j++)
      if (i == 1)
        if (j == grX.size / 2)
          printf("     1    | %4g |      --      |  %2.9lf  |      --      \n", grX.nodes[j], func(grX.nodes[j]));
        else
          printf("          | %4g |      --      |  %2.9lf  |      --      \n", grX.nodes[j], func(grX.nodes[j]));
      else
      {
        if (j == grX.size / 2)
          printf("%6i    | %4g | %2.9lf  |  %2.9lf  | %g%%\n", i, grX.nodes[j], a = grY2.nodes[j], b = grY.nodes[j], fabs(a - b) / b * 100);

        printf("          | %4g | %2.9lf  |  %2.9lf  | %g%%\n", grX.nodes[j], a = grY2.nodes[j], b = grY.nodes[j], fabs(a - b) / b * 100);
      }

    printf("-------------------------------------------------------------\n");
    GridFree(&grY);
  }

  GridFree(&grX);
#endif
  /* Predictor-corrector */
#if 0
  printf("solution by predictor-corrector method (4 order) for d.e.: x2y' - 2xy = 3y, x = [1, 4]\nintervals - num of intervals between nodes\n\n");
  printf("intervals | node |   Pred-corr  |    exact    | error\n");
  printf("-------------------------------------------------------------\n");
  for (i = 1; i <= 50; i += 7)
  {
    double a, b;
    grY = AdamsPredCorr(diffEq, func(grX.nodes[0]), &grX, i);

    for (j = 0; j < grX.size; j++)
      if (i == 1)
        if (j == grX.size / 2)
          printf("     1    | %4g |      --      | %2.9lf |      --      \n", grX.nodes[j], func(grX.nodes[j]));
        else
          printf("          | %4g |      --      | %2.9lf |      --      \n", grX.nodes[j], func(grX.nodes[j]));
      else
      {
        if (j == grX.size / 2)
          printf("%6i    | %4g | %2.9lf  | %2.9lf | %g\n", i, grX.nodes[j], a = grY.nodes[j], b = func(grX.nodes[j]), fabs(a - b));

        printf("          | %4g | %2.9lf  | %2.9lf | %g\n", grX.nodes[j], a = grY.nodes[j], b = func(grX.nodes[j]), fabs(a - b));
      }

    printf("-------------------------------------------------------------\n");
    GridFree(&grY);
  }

  printf("\n\nerror in start condition = 1%%\n\n");
  printf("intervals | node |   Pred-corr  | without error | difference (%%)\n");
  printf("-------------------------------------------------------------\n");
  for (i = 1; i <= 50; i += 7)
  {
    double a, b;
    grY = AdamsPredCorr(diffEq, func(grX.nodes[0]), &grX, i);
    grY2 = AdamsPredCorr(diffEq, (a = func(grX.nodes[0])) + 0.01 * a, &grX, i);

    for (j = 0; j < grX.size; j++)
      if (i == 1)
        if (j == grX.size / 2)
          printf("     1    | %4g |      --      |  %2.9lf  |      --      \n", grX.nodes[j], func(grX.nodes[j]));
        else
          printf("          | %4g |      --      |  %2.9lf  |      --      \n", grX.nodes[j], func(grX.nodes[j]));
      else
      {
        if (j == grX.size / 2)
          printf("%6i    | %4g | %2.9lf  |  %2.9lf  | %g%%\n", i, grX.nodes[j], a = grY2.nodes[j], b = grY.nodes[j], fabs(a - b) / b * 100);

        printf("          | %4g | %2.9lf  |  %2.9lf  | %g%%\n", grX.nodes[j], a = grY2.nodes[j], b = grY.nodes[j], fabs(a - b) / b * 100);
      }

    printf("-------------------------------------------------------------\n");
    GridFree(&grY);
  }

  GridFree(&grX);
#endif
  getchar();
}
