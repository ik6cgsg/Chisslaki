#include "rootf.h"

int HalfDiv( double (*f)(double), double left, double right, double *x, double TRESHOLD )
{
  double tmp = f(2);
  int iter_num = 0;

  if (f(left) * f(right) > 0)
    return iter_num;

  while (right - left > TRESHOLD)
  {
    tmp = (right + left) / 2;
    if (f(tmp) * f(right) <= 0)
      left = tmp;
    else if (f(left) * f(tmp) <= 0)
      right = tmp;
    iter_num++;
  }

  *x = (right + left) / 2;

  return iter_num;
}

int SimpleIter( double (*f)(double), double x0, double q, double *x, double TRESHOLD )
{
  double x_k = f(x0), x_k_1 = x0;
  int iter = 1;

  while (f(x_k) != 0 && fabs(x_k - x_k_1) > (1 - q) / q * TRESHOLD)
  {
    x_k_1 = x_k;
    x_k = f(x_k);
    iter++;
  }

  *x = x_k;

  return iter;
}
