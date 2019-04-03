#include "gradient.h"

vec Grad( matr * a, vec * b, vec * x )
{
  vec ax = MatrMulVec(a, x), ax_b, grad;

  ax_b = VecMinusVec(&ax, b);
  grad = VecMulNum(&ax_b, 2);

  FreeVec(&ax);
  FreeVec(&ax_b);

  return grad;
}

static double s_getAlpha( vec *gk, matr *a )
{
  vec ag = MatrMulVec(a, gk);
  double al = Scalar(gk, gk) / (2 * Scalar(&ag, gk));

  FreeVec(&ag);

  return al;
}

vec SolveGrad( matr * a, vec * b, vec * x0, double eps, int *iter )
{
  vec xk_1 = Zeros(a->size), xk = *x0, dif = VecMinusVec(&xk, &xk_1);
  int start = 1;
  *iter = 0;

  FreeVec(&xk_1);

  while (NormEucl(&dif) > eps)
  {
    vec gk = Grad(a, b, &xk), agk;
    double alpha = s_getAlpha(&gk, a);

    FreeVec(&dif);

    agk = VecMulNum(&gk, alpha);
    xk_1 = VecMinusVec(&xk, &agk);

    FreeVec(&agk);
    FreeVec(&gk);

    dif = VecMinusVec(&xk, &xk_1);

    if (!start)
      FreeVec(&xk);
    else
      start = 0;

    xk = xk_1;

    *iter += 1;
  }

  FreeVec(&dif);

  return xk;
}
