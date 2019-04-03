#include "interpol.h"

double static s_li( grid const *gr, uint i, double x )
{
  double prod = 1;
  uint j;

  for (j = 0; j < gr->size; j++)
    if (j != i)
      prod *= (x - gr->nodes[j]) / (gr->nodes[i] - gr->nodes[j]);

  return prod;
} /* End of 's_li' function */

double InterPolLagrange( double (*f)(double), grid const *gr, double x )
{
  double sum = 0;
  uint i;

  for (i = 0; i < gr->size; i++)
    sum += f(gr->nodes[i]) * s_li(gr, i, x);

  return sum;
} /* End of 'InterPolLagrange' function */

void GenerateY( grid const *x, grid *y, double (*f)(double), uint nodes_num, grid_type gt )
{
  uint i;
  grid inter_gr;

  GridInit(&inter_gr, nodes_num);

  switch (gt)
  {
  case GRID_UNIFORM:
    GridUniform(&inter_gr, x->nodes[0], x->nodes[x->size - 1]);
    break;
  case GRID_CHEBYSHEV:
    GridChebyshev(&inter_gr, x->nodes[0], x->nodes[x->size - 1]);
    break;
  case GRID_USER:
    GridHilb(&inter_gr, x->nodes[0], x->nodes[x->size - 1]);
    break;
  default:
    GridUniform(&inter_gr, x->nodes[0], x->nodes[x->size - 1]);
    break;
  }

  for (i = 0; i < x->size; i++)
    y->nodes[i] = InterPolLagrange(f, &inter_gr, x->nodes[i]);

  GridFree(&inter_gr);
} /* End of 'GenerateY' function */
