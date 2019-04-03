#include "eigen.h"
#include "rootf.h"
#include "decompos.h"

struct
{
  matr m;
} Curm;

static double s_getProd(int i, int size)
{
  int j;
  double prod = 1;

  for (j = 0; j <= i; j++)
    prod *= Curm.m.a[size - j - 1][size - j - 2];

  return prod;
}

static double s_recCharactHess(double x, int size)
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

static double s_charactDetHess(double x)
{
  return s_recCharactHess(x, Curm.m.size);
}

vec FindHessEigenValues(matr *m, double step, double eps, int *iter)
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

vec InverseIter(matr *m, vec *y0, double eps, int *iter)
{
  vec
    eigs = Zeros(2),
    y_i = Zeros(y0->length),
    y_i_n = Zeros(y0->length),
    y_i_n_1 = Zeros(y0->length),
    y_k, gw, y_k_n, dif = Ones(m->size);
  double mu_i, mu_k = 0, mu_i_1, gamma;
  matr q, r;

  *iter = 0;

  Curm.m = Unit(m->size);
  CopyMatr(&Curm.m, m);

  if (!QRDecompose(&Curm.m, &q, &r))
  {
    FreeVec(&y_i);
    FreeVec(&y_i_n);
    return eigs;
  }

  /* Start conditions */
  CopyVec(&y_i, y0);
  mu_i = NormEucl(y0);
  mu_i_1 = mu_i + 2 * eps;    // just to step into loop

  /* ^ norm. method for first dominant eigenvalue */
  while (/*fabs(mu_i - mu_i_1)*/NormEucl(&dif) > eps)
  {
    CopyVec(&y_i_n_1, &y_i_n);
    FreeVec(&y_i_n);
    y_i_n = VecMulNum(&y_i, 1 / mu_i);
    FreeVec(&y_i);
    y_i = Solve(&Curm.m, &y_i_n, &q, &r);

    FreeVec(&dif);
    dif = VecMinusVec(&y_i_n, &y_i_n_1);

    mu_i_1 = mu_i;
    mu_i = NormEucl(&y_i);

    (*iter)++;
  }

  /////// y_i - close to last eigenvector
  CopyVec(&y_i_n, &y_i);

  /////// ortogonalization
  gamma = Scalar(y0, &y_i_n);
  gw = VecMulNum(&y_i_n, gamma);
  y_k = VecMinusVec(y0, &gw);
  y_k_n = Zeros(y_k.length);

  /* Start conditions */
  mu_k = NormEucl(&y_k);
  mu_i_1 = mu_k + 2 * eps;    // just to step into loop
  FreeVec(&dif);
  dif = Ones(m->size);

  /* ^ norm. + ort. method for second dominant eigenvalue */
  while (/*fabs(mu_k - mu_i_1)*/NormEucl(&dif) > eps)
  {
    /* ortogonalization */
    if ((*iter) % 2 == 0)
    {
      gamma = Scalar(&y_k, &y_i_n);
      FreeVec(&gw);
      gw = VecMulNum(&y_i_n, gamma);
      FreeVec(&y_i);
      y_i = VecMinusVec(&y_k, &gw);
      CopyVec(&y_k, &y_i);
    }

    /* ^ norm. method */
    CopyVec(&y_i_n_1, &y_k_n);
    FreeVec(&y_k_n);
    y_k_n = VecMulNum(&y_k, 1 / mu_k);
    FreeVec(&y_k);
    y_k = Solve(&Curm.m, &y_k_n, &q, &r);

    FreeVec(&dif);
    dif = VecMinusVec(&y_k_n, &y_i_n_1);

    mu_i_1 = mu_k;
    mu_k = NormEucl(&y_k);

    (*iter)++;
  }

  /* First & second minimal eigenvalues */
  eigs.v[0] = 1 / mu_i;
  eigs.v[1] = 1 / mu_k;

  FreeMatr(&Curm.m);
  FreeMatr(&q);
  FreeMatr(&r);
  FreeVec(&y_i);
  FreeVec(&y_i_n);
  FreeVec(&y_k_n);
  FreeVec(&y_k);
  FreeVec(&gw);

  return eigs;
}
