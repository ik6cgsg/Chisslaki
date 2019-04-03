#include <math.h>

#include "diff.h"

#define EPSYLON 1e-9

static int equal(double a, double b)
{
  if (fabs(a - b) < EPSYLON)
    return 1;
  return 0;
}

grid RungeKutta(funcptr2 diffEq, double y0, grid * gr, uint numInt)
{
  grid res, local;
  double hGr = gr->nodes[1] - gr->nodes[0];
  double h = hGr / numInt;
  uint i, j;

  GridInit(&res, gr->size);
  res.nodes[0] = y0;

  GridInit(&local, numInt + 1);

  for (i = 1; i < gr->size; i++)
  {
    double y_i = res.nodes[i - 1];
    double eta1, eta2, eta3, eta4, dy;
    GridUniform(&local, gr->nodes[i - 1], gr->nodes[i]);

    for(j = 0; j < local.size - 1; j++)
    {
      eta1 = diffEq(local.nodes[j], y_i);
      eta2 = diffEq(local.nodes[j] + h / 2, y_i + h / 2 * eta1);
      eta3 = diffEq(local.nodes[j] + h / 2, y_i + h / 2 * eta2);
      eta4 = diffEq(local.nodes[j] + h, y_i + h * eta3);
      dy = h / 6 * (eta1 + 2 * eta2 + 2 * eta3 + eta4);
      y_i += dy;
    }

    res.nodes[i] = y_i;
  }

  GridFree(&local);
  return res;
}

grid AdamsBashforth(funcptr2 diffEq, double y0, grid * gr, uint numInt)
{
  grid res, tmpY, tmpX;
  double hGr = gr->nodes[1] - gr->nodes[0];
  double h = hGr / numInt;
  double y_j, y_j_1, y_j_2, y_j_3;
  double x_j_3 = gr->nodes[0], x_j_2 = x_j_3 + h, x_j_1 = x_j_2 + h, x_j = x_j_1 + h;
  uint i, j;

  GridInit(&res, gr->size);
  res.nodes[0] = y0;

  GridInit(&tmpX, 4);
  tmpX.nodes[0] = x_j_3;
  tmpX.nodes[1] = x_j_2;
  tmpX.nodes[2] = x_j_1;
  tmpX.nodes[3] = x_j;

  tmpY = RungeKutta(diffEq, y0, &tmpX, 1);
  y_j_3 = y0;
  y_j_2 = tmpY.nodes[1];
  y_j_1 = tmpY.nodes[2];
  y_j = tmpY.nodes[3];

  GridFree(&tmpX);
  GridFree(&tmpY);

  for (i = 1; i < gr->size; i++)
  {
    double y_i = res.nodes[i - 1];

    if (i == 1)
      j = 3;
    else 
      j = 0;

    for(j; j < numInt; j++)
    {
      y_i = y_j + h / 24 * (55 * diffEq(x_j, y_j) - 59 * diffEq(x_j_1, y_j_1) +
                            37 * diffEq(x_j_2, y_j_2) - 9 * diffEq(x_j_3, y_j_3));
      y_j_3 = y_j_2;
      y_j_2 = y_j_1;
      y_j_1 = y_j;
      y_j = y_i;

      x_j_3 = x_j_2;
      x_j_2 = x_j_1;
      x_j_1 = x_j;
      x_j += h;
    }

    res.nodes[i] = y_i;
  }

  return res;
}

grid AdamsPredCorr(funcptr2 diffEq, double y0, grid * gr, uint numInt)
{
  grid res, tmpX, tmpY;
  double hGr = gr->nodes[1] - gr->nodes[0];
  double h = hGr / numInt;
  double y_j, y_j_1, y_j_2, y_j_3;
  double x_j_3 = gr->nodes[0], x_j_2 = x_j_3 + h, x_j_1 = x_j_2 + h, x_j = x_j_1 + h;
  uint i, j;

  GridInit(&res, gr->size);
  res.nodes[0] = y0;

  GridInit(&tmpX, 4);
  tmpX.nodes[0] = x_j_3;
  tmpX.nodes[1] = x_j_2;
  tmpX.nodes[2] = x_j_1;
  tmpX.nodes[3] = x_j;

  tmpY = RungeKutta(diffEq, y0, &tmpX, 1);
  y_j_3 = y0;
  y_j_2 = tmpY.nodes[1];
  y_j_1 = tmpY.nodes[2];
  y_j = tmpY.nodes[3];

  GridFree(&tmpX);
  GridFree(&tmpY);

  for (i = 1; i < gr->size; i++)
  {
    double y_i = res.nodes[i - 1];

    if (i == 1)
      j = 3;
    else 
      j = 0;

    for (j; j < numInt; j++)
    {
      y_i = y_j + h / 24 * (55 * diffEq(x_j, y_j) - 59 * diffEq(x_j_1, y_j_1) +
                            37 * diffEq(x_j_2, y_j_2) - 9 * diffEq(x_j_3, y_j_3));

      y_i = y_j + h / 24 * (9 * diffEq(x_j + h, y_i) + 19 * diffEq(x_j, y_j) -
                            5 * diffEq(x_j_1, y_j_1) + diffEq(x_j_2, y_j_2));
      y_j_3 = y_j_2;
      y_j_2 = y_j_1;
      y_j_1 = y_j;
      y_j = y_i;

      x_j_3 = x_j_2;
      x_j_2 = x_j_1;
      x_j_1 = x_j;
      x_j += h;
    }

    res.nodes[i] = y_i;
  }
  return res;
}

struct
{
  funcptr p, q, f;
  grid *delta, *gamma, *deltax, *gammax;
} Ptrs;

static double deltaFunc( double x )
{
  uint start = 0, end = Ptrs.deltax->size;
  uint ind = (start + end) >> 1;

  while (!equal(Ptrs.deltax->nodes[ind], x) && ind >= 0 && ind < Ptrs.deltax->size)
  {
    if (Ptrs.deltax->nodes[ind] > x)
      end = ind;
    else
      start = ind;

    ind = (start + end) >> 1;
  }

  if (ind < 0 || ind >= Ptrs.deltax->size)
    return 0;

  return Ptrs.delta->nodes[ind];
}

static double gammaFunc( double x )
{
  uint start = 0, end = Ptrs.gammax->size;
  uint ind = (start + end) >> 1;

  while (!equal(Ptrs.gammax->nodes[ind], x) && ind >= 0 && ind < Ptrs.gammax->size)
  {
    if (Ptrs.gammax->nodes[ind] > x)
      end = ind;
    else
      start = ind;

    ind = (start + end) >> 1;
  }

  if (ind < 0 || ind >= Ptrs.gammax->size)
    return 0;

  return Ptrs.gamma->nodes[ind];
}

static double deltaEq( double x, double delta )
{
  return -Ptrs.p(x) * delta - delta * delta - Ptrs.q(x);
}

static double gammaEq( double x, double gamma )
{
  return -(Ptrs.p(x) + deltaFunc(x)) * gamma + Ptrs.f(x);
}

static double diff( double x, double y )
{
  return deltaFunc(x) * y + gammaFunc(x);
}

grid DiffRun(funcptr p, funcptr q, funcptr f, grid * gr, double stCond[6], uint numInt)
{
  grid res, delta, gamma, deltax, gammax;
  GridInit(&res, gr->size);
  double alpha0 = stCond[0], alpha1 = stCond[1], A = stCond[2];
  double beta0 = stCond[3], beta1 = stCond[4], B = stCond[5];
  double b = gr->nodes[gr->size - 1];

  Ptrs.p = p;
  Ptrs.q = q;
  Ptrs.f = f;

  GridInit(&deltax, (gr->size - 1) * numInt * 2 * numInt * 2 + 1);
  GridUniform(&deltax, gr->nodes[0], gr->nodes[gr->size - 1]);
  Ptrs.deltax = &deltax;

  delta = RungeKutta(deltaEq, -alpha0 / alpha1, &deltax, numInt);
  Ptrs.delta = &delta;

  GridInit(&gammax, (gr->size - 1) * numInt * 2 + 1);
  GridUniform(&gammax, gr->nodes[0], gr->nodes[gr->size - 1]);
  Ptrs.gammax = &gammax;

  gamma = RungeKutta(gammaEq, A / alpha1, &gammax, numInt);
  Ptrs.gamma = &gamma;

  GridReverse(gr);
  res = RungeKutta(diff, (B - beta1 * gammaFunc(b)) / (beta0 + beta1 * deltaFunc(b)), gr, numInt);
  GridReverse(gr);

  GridReverse(&res);

  GridFree(&delta);
  GridFree(&gamma);
  GridFree(&deltax);
  GridFree(&gammax);

  return res;
}
