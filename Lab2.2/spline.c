#include <stdlib.h>

#include "spline.h"

void SplineInit( spline * sp, uint size, double begin, double end )
{
  GridInit(&(sp->gr), size);
  GridUniform(&(sp->gr), begin, end);

  sp->polys = malloc((size - 1) * sizeof(cub_poly));
} /* End of 'SplineInit' function */

void BuildNaturalSpline( spline * sp, double(*f)(double) )
{
  uint i, n = sp->gr.size - 1;
  double *delta, *lam;
  double h = sp->gr.nodes[1] - sp->gr.nodes[0];

  delta = malloc(sizeof(double) * (n - 1));
  lam = malloc(sizeof(double) * (n - 1));

  delta[0] = -1.0 / 4;
  lam[0] = 3 / (4 * h * h) * (f(sp->gr.nodes[2]) - 2 * f(sp->gr.nodes[1]) + f(sp->gr.nodes[0]));

  for (i = 2; i < n; i ++)
  {
    delta[i - 1] = -1 / (4 + delta[i - 2]);
    lam[i - 1] = (3 / (h * h) * (f(sp->gr.nodes[i + 1]) - 2 * f(sp->gr.nodes[i]) + f(sp->gr.nodes[i - 1])) - lam[i - 2]) / (4 + delta[i - 2]);
  }

  sp->polys[n - 1].c = 0;

  for (i = n - 1; i > 0; i--)
    sp->polys[i - 1].c = delta[i - 1] * sp->polys[i].c + lam[i - 1];

  for (i = 0; i < n; i++)
  {
    double ck_1;

    if (i == 0)
      ck_1 = 0;
    else
      ck_1 = sp->polys[i - 1].c;

    sp->polys[i].b = (f(sp->gr.nodes[i + 1]) - f(sp->gr.nodes[i])) / h + h / 3 * (2 * sp->polys[i].c + ck_1);
    sp->polys[i].d = (sp->polys[i].c - ck_1) / 3 / h;
    sp->polys[i].a = f(sp->gr.nodes[i + 1]);
  }

  free(delta);
  free(lam);
} /* End of 'BuildNaturalSpline' function */

void BuildSecDerivSpline( spline * sp, double(*f)(double), double sec_deriv_x0, double sec_deriv_xn )
{
  uint i, n = sp->gr.size - 1;
  double *delta, *lam;
  double h = sp->gr.nodes[1] - sp->gr.nodes[0];

  delta = malloc(sizeof(double) * (n - 1));
  lam = malloc(sizeof(double) * (n - 1));

  delta[0] = -1.0 / 4;
  lam[0] = 3 / (4 * h * h) * (f(sp->gr.nodes[2]) - 2 * f(sp->gr.nodes[1]) + f(sp->gr.nodes[0]));

  for (i = 2; i < n; i ++)
  {
    delta[i - 1] = -1 / (4 + delta[i - 2]);
    lam[i - 1] = (3 / (h * h) * (f(sp->gr.nodes[i + 1]) - 2 * f(sp->gr.nodes[i]) + f(sp->gr.nodes[i - 1])) - lam[i - 2]) / (4 + delta[i - 2]);
  }

  sp->polys[n - 1].c = sec_deriv_xn / 2;

  for (i = n - 1; i > 0; i--)
    sp->polys[i - 1].c = delta[i - 1] * sp->polys[i].c + lam[i - 1];

  for (i = 0; i < n; i++)
  {
    double ck_1;

    if (i == 0)
      ck_1 = sec_deriv_x0 / 2;
    else
      ck_1 = sp->polys[i - 1].c;

    sp->polys[i].b = (f(sp->gr.nodes[i + 1]) - f(sp->gr.nodes[i])) / h + h / 3 * (2 * sp->polys[i].c + ck_1);
    sp->polys[i].d = (sp->polys[i].c - ck_1) / 3 / h;
    sp->polys[i].a = f(sp->gr.nodes[i + 1]);
  }

  free(delta);
  free(lam);
} /* End of 'BuildSecDerivSpline' function */

double GetValue( spline * sp, double x )
{
  uint a = 0, b = sp->gr.size - 1;

  while (a < b)
  {
    uint mid = a + (b - a) / 2;

    if (x > sp->gr.nodes[mid])
      a = mid + 1;
    else
      b = mid;
  }

  if (a == 0)
    a++;

  return sp->polys[a - 1].a + sp->polys[a - 1].b * (x - sp->gr.nodes[a]) +
         sp->polys[a - 1].c * pow((x - sp->gr.nodes[a]), 2) +
         sp->polys[a - 1].d * pow((x - sp->gr.nodes[a]), 3);
} /* End of 'GetValue' function */

void SplineFree( spline * sp )
{
  GridFree(&(sp->gr));
  free(sp->polys);
} /* End of 'SplineFree' function */
