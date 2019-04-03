#include <stdlib.h>

#include "grid.h"

#define PI 3.14159265358979323

int GridInit( grid *gr, uint size )
{
  if ((gr->nodes = malloc(size * sizeof(double))) == NULL)
    return 0;

  gr->size = size;

  return 1;
} /* End of 'GridInit' function */

void GridUniform( grid *gr, double begin, double end )
{
  double h = (end - begin) / (gr->size - 1);
  uint i;

  for (i = 0; i < gr->size; i++)
    gr->nodes[i] = begin + i * h;
} /* End of 'GridUniform' function */

void GridChebyshev( grid *gr, double begin, double end )
{
  uint i;

  for (i = 0; i < gr->size; i++)
    gr->nodes[gr->size - i - 1] = (begin + end) / 2 + (end - begin) / 2 * cos(1.0 * (2 * (i + 1) - 1) / (2 * gr->size) * PI);
} /* End of 'GridChebyshev' function */

void GridHilb( grid *gr, double begin, double end )
{
  uint i;

  gr->nodes[0] = begin;
  gr->nodes[gr->size - 1] = end;

  for (i = gr->size - 2; i >= gr->size / 2; i--)
  {
    gr->nodes[gr->size - i - 1] = begin + 1.0 / (i + 1) * (end - begin);
    gr->nodes[i] = end - 1.0 / (i + 1) * (end - begin);
  }
} /* End of 'GridRandom' function */

void GridFree( grid *gr )
{
  if (gr->size != 0)
    free(gr->nodes);

  gr->size = 0;
} /* End of 'GridFree' function */
