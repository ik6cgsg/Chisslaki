#include <math.h>

#include "diff.h"

#define EPSYLON 1e-7

static int equal(double a, double b)
{
  if (fabs(a - b) < EPSYLON)
    return 1;
  return 0;
}

grid RungeKutta(double(*diffEq)(double, double), double y0, grid * gr, uint numInt)
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

grid AdamsBashforth(double(*diffEq)(double, double), double y0, grid * gr, uint numInt)
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

grid AdamsPredCorr(double(*diffEq)(double, double), double y0, grid * gr, uint numInt)
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
