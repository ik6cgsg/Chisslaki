#include <stdio.h>
#include <math.h>

#include "diff.h"

#define e 2.718281828459045
#define a 1
#define b 1.25
#define size_t 4

double p(double x)
{
  return -(2 * x + 1) / x;
}

double q(double x)
{
  return (x + 2) / x;
}

double f(double x)
{
  return x * exp(x);
}

double y(double x)
{
  return x * x * exp(x);
}

double dy(double x)
{
  return exp(x) * (2 * x + x * x);
}

void main(void)
{
  grid gx, gy;
  //double al0 = 100, al1 = 10000, bt0 = 10000, bt1 = 0;
  double al0 = 1, al1 = 1, bt0 = 1, bt1 = 1;
  double cond[] = {al0, al1, al0 * y(a) + al1 * dy(a), bt0, bt1, bt0 * y(b) + bt1 * dy(b)};
  int i;

  GridInit(&gx, size_t);
  GridUniform(&gx, a, b);

  printf("solution by Differential run for d.e.: xy'' - (2x + 1)y' + (x + 2)y = x^2e^x, x = [1, 1.25]\nintervals - num of intervals between nodes\n");
  printf("start conditions: a0 = %lf, a1 = %lf, A = %lf\nb0 = %lf, b1 = %lf, B = %lf\n\n", cond[0], cond[1], cond[2], cond[3], cond[4], cond[5]);
  printf("intervals |   node  |  Diff run   |    exact    | error\n");
  printf("-------------------------------------------------------------\n");
  for (i = 1; i <= 50; i += 7)
  {
    double c, d;
    int j;

    gy = DiffRun(p, q, f, &gx, cond, i);

    for (j = 0; j < gx.size; j++)
    {
      if (j == gx.size / 2)
        printf("%6i    | %2.5lf | %2.9lf | %2.9lf | %g\n", i, gx.nodes[j], c = gy.nodes[j], d = y(gx.nodes[j]), fabs(c - d));
      else
        printf("          | %2.5lf | %2.9lf | %2.9lf | %g\n", gx.nodes[j], c = gy.nodes[j], d = y(gx.nodes[j]), fabs(c - d));
    }

    printf("-------------------------------------------------------------\n");
    GridFree(&gy);
  }


  printf("\n\nerror in start condition = 1%%\n\n");
  printf("intervals |  node   |   Diff run   | without error | difference (%%)\n");
  printf("----------------------------------------------------------------\n");
  for (i = 1; i <= 50; i += 7)
  {
    double c, d;
    int j;
    double cond2[6];
    grid gy2;

    for (j = 0; j < 6; j++)
      cond2[j] = cond[j];

    cond2[2] *= 1.01;
    cond2[5] *= 1.01;

    gy = DiffRun(p, q, f, &gx, cond, i);
    gy2 = DiffRun(p, q, f, &gx, cond2, i);

    for (j = 0; j < gx.size; j++)
    {
      if (j == gx.size / 2)
        printf("%6i    | %2.5lf | %2.9lf  |  %2.9lf  | %g%%\n", i, gx.nodes[j], c = gy2.nodes[j], d = gy.nodes[j], fabs(c - d) / d * 100);
      else
        printf("          | %2.5lf | %2.9lf  |  %2.9lf  | %g%%\n", gx.nodes[j], c = gy2.nodes[j], d = gy.nodes[j], fabs(c - d) / d * 100);
    }

    printf("----------------------------------------------------------------\n");
    GridFree(&gy);
    GridFree(&gy2);
  }

  GridFree(&gx);

  getchar();
}
