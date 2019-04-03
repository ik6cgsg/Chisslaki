#include "integr.h"
#include "grid.h"

static double st_integrate( double(*f)(double), grid *gr )
{
  return 3.0 / 8 * (gr->nodes[1] - gr->nodes[0]) * (f(gr->nodes[0]) + 3 * f(gr->nodes[1]) + 3 * f(gr->nodes[2]) + f(gr->nodes[3]));
} /* End of 'st_integrate' function */

double Integrate( double(*f)(double), double a, double b, int pNum )
{
  int i;
  grid gr, local;
  double sum = 0;

  GridInit(&gr, pNum);
  GridUniform(&gr, a, b);

  GridInit(&local, 4);

  for (i = 0; i < pNum - 1; i++)
  {
    GridUniform(&local, gr.nodes[i], gr.nodes[i + 1]);

    sum += st_integrate(f, &local);
  }

  GridFree(&gr);
  GridFree(&local);
  
  return sum;
} /* End of 'Integrate' function */

double IntegrateRunge( double(*f)(double), double a, double b, double eps, int *pNum )
{
  int N = 1;
  double s1, s2;

  do
  {
    N++;
    s1 = Integrate(f, a, b, N);
    s2 = Integrate(f, a, b, 2 * N);
  } while (fabs(s2 - s1) / 80 > eps);

  *pNum = 2 * N;
  return s2;
} /* End of 'IntegrateRunge' function */

double NewLeyb( double(*fInt)(double), double a, double b )
{
  return fInt(b) - fInt(a);
} /* End of 'NewLeyb' function */
