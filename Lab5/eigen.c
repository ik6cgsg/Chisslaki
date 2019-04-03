#include "eigen.h"
#include "rootf.h"

struct
{
  matr m;
} Curm;

static double s_getProd( int i, int size )
{
  int j;
  double prod = 1;

  for (j = 0; j <= i; j++)
    prod *= Curm.m.a[size - j - 1][size - j - 2];

  return prod;
}

static double s_recCharactHess( double x, int size )
{
  double res = 0;
  int i;

  if (size == 0)
    return 1;

  res += (Curm.m.a[size - 1][size - 1] - x) * s_recCharactHess(x, size - 1);

  for (i = 0; i < size - 1; i++)
    res += pow(-1, i + 1) * Curm.m.a[size - i - 2][size - 1] * s_getProd(i, size) * s_recCharactHess(x, size - 2 - i);

  return res;
}

static double s_charactDetHess( double x )
{
  return s_recCharactHess(x, Curm.m.size);
}

vec FindHessEigenValues( matr *m, double step, double eps, int *iter )
{
  double res, x_start = -InfNorm(m), x_end = InfNorm(m), x_k = x_start, x_k_1 = x_k;
  int i = 0;
  vec eigens = Zeros(m->size);
  *iter = 0;

  Curm.m = ToHess(m);

  while (x_k < x_end - step)
  {
    x_k_1 = x_k;
    x_k += step;
    if (s_charactDetHess(x_k_1) * s_charactDetHess(x_k) < 0)
    {
      HalfDiv(s_charactDetHess, x_k_1, x_k, &res, eps);
      eigens.v[i++] = res;
    }
    (*iter)++;
  }

  FreeMatr(&Curm.m);

  return eigens;
}
