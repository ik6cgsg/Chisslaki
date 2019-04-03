#include "integr.h"
#include "grid.h"

static double st_integrate38( double(*f)(double), grid *gr )
{
  return 3.0 / 8 * (gr->nodes[1] - gr->nodes[0]) * (f(gr->nodes[0]) + 3 * f(gr->nodes[1]) + 3 * f(gr->nodes[2]) + f(gr->nodes[3]));
} /* End of 'st_integrate' function */

static double st_integrate_gauss( double(*f)(double), double a, double b )
{
  double A[4] = {0.347855, 0.652145, 0.652145, 0.347855};
  double t[4] = {-0.861136, -0.339981, 0.339981, 0.861136};
  int i;
  double sum = 0;

  for (i = 0; i < 4; i++)
    sum += A[i] * f((a + b) / 2 + (b - a) / 2 * t[i]);

  return sum * (b - a) / 2;
} /* End of 'st_integrate' function */

double Integrate38( double(*f)(double), double a, double b, int pNum )
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

    sum += st_integrate38(f, &local);
  }

  GridFree(&gr);
  GridFree(&local);
  
  return sum;
} /* End of 'Integrate' function */

double IntegrateGauss(double(*f)(double), double a, double b, int pNum)
{
  int i;
  grid gr;
  double sum = 0;

  GridInit(&gr, pNum);
  GridUniform(&gr, a, b);

  for (i = 0; i < pNum - 1; i++)
    sum += st_integrate_gauss(f, gr.nodes[i], gr.nodes[i + 1]);

  GridFree(&gr);

  return sum;
} /* End of 'IntegrateGauss' function */

double Integrate38ByEps(double(*f)(double), double a, double b, double eps, int * pNum)
{
  int N = 1;
  double s1, s2;

  do
  {
    N++;
    s1 = Integrate38(f, a, b, N);
    s2 = Integrate38(f, a, b, 2 * N);
  } while (fabs(s2 - s1) / 15 > eps);

  *pNum = 2 * N;
  return s2;
} /* End of 'Integrate38ByEps' function */

double IntegrateGaussByEps(double(*f)(double), double a, double b, double eps, double res, int * pNum)
{
  int N = 1;
  double s;

  do
  {
    N *= 2;
    s = IntegrateGauss(f, a, b, N);
  } while (fabs(res - s) > eps);

  *pNum = N;
  return s;
} /* End of 'IntegrateGaussByEps' function */

double NewLeib( double(*fInt)(double), double a, double b )
{
  return fInt(b) - fInt(a);
} /* End of 'NewLeyb' function */
