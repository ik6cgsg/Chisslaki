#include "matr.h"
#include "vec.h"

#define ZERO(A, EPS) (A > -EPS && A < EPS)

static double s_determ2x2(matr *m)
{
  return m->a[0][0] * m->a[1][1] - m->a[0][1] * m->a[1][0];
}

static int s_triangle(matr *m)
{
  int i, j;

  for (i = 0; i < m->size; i++)
    for (j = 0; j < i; j++)
      if (m->a[i][j] != 0)
        return 0;
  return 1;
}

void SwapPtr(double **a, double **b)
{
  double *tmp = *a;
  *a = *b;
  *b = tmp;
}

static matr s_alg(matr *m, int x, int y)
{
  int i, j;
  matr m_new = Unit(m->size);

  for (i = 0; i < m->size; i++)
    for (j = 0; j < m->size; j++)
      m_new.a[i][j] = i == y && j == x ? 1 : i == y || j == x ? 0 : m->a[i][j];

  return m_new;
}

static void s_multNum(matr *m, double n)
{
  int i, j;

  for (i = 0; i < m->size; i++)
    for (j = 0; j < m->size; j++)
      m->a[i][j] *= n;
}

matr Unit(int size)
{
  matr m;
  int i, j;

  m.size = size;
  m.a = malloc(sizeof(double *) * size);

  for (i = 0; i < size; i++)
    m.a[i] = malloc(sizeof(double) * size);

  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      m.a[i][j] = i == j ? 1 : 0;

  return m;
}

matr RandMatr(int size)
{
  matr m;
  int i, j;

  m.size = size;
  m.a = malloc(sizeof(double *) * size);

  for (i = 0; i < size; i++)
    m.a[i] = malloc(sizeof(double) * size);

  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      m.a[i][j] = (int)RAND_VALUE(MAX_RAND_VALUE);

  return m;
}

matr Hilb(int size)
{
  matr m = Unit(size);
  int i, j;

  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      m.a[i][j] = 1.0 / (i + j + 1);

  return m;
}

matr MakeMatr(int size, double *a)
{
  int i, j;
  matr mn = Unit(size);

  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      mn.a[i][j] = a[i * size + j];

  return mn;
}

matr MatrMulNum(matr * m, double n)
{
  int i, j;
  matr mn = Unit(m->size);

  for (i = 0; i < m->size; i++)
    for (j = 0; j < m->size; j++)
      mn.a[i][j] = m->a[i][j] * n;

  return mn;
}

matr MatrMulMatr(matr *m1, matr *m2, int ind)
{
  matr prod = Unit(m1->size);
  int i, j;

  for (i = ind; i < prod.size; i++)
    for (j = ind; j < prod.size; j++)
    {
      double sum = 0;
      int k;

      for (k = ind; k < prod.size; k++)
        sum += m1->a[i][k] * m2->a[k][j];

      prod.a[i][j] = sum;
    }

  return prod;
}

matr MatrPlusMatr(matr * m1, matr * m2)
{
  int i, j;
  matr m = Unit(m1->size);

  for (i = 0; i < m1->size; i++)
    for (j = 0; j < m1->size; j++)
      m.a[i][j] = m1->a[i][j] + m2->a[i][j];

  return m;
}

matr MatrMinusMatr(matr * m1, matr * m2)
{
  int i, j;
  matr m = Unit(m1->size);

  for (i = 0; i < m1->size; i++)
    for (j = 0; j < m1->size; j++)
      m.a[i][j] = m1->a[i][j] - m2->a[i][j];

  return m;
}

matr InverseMatr(matr * m)
{
  int i, j;
  double det;
  matr m_t = Unit(m->size);

  for (i = 0; i < m->size; i++)
    for (j = 0; j < m->size; j++)
    {
      /* Building algebraic transpose adjunct */
      matr tmp = s_alg(m, i, j);

      m_t.a[i][j] = Det(&tmp);
      FreeMatr(&tmp);
    }

  det = Det(m);

  s_multNum(&m_t, 1 / det);

  return m_t;
}

matr OffsetMatr(matr * m)
{
  matr m_off = Unit(m->size);
  int i, j;

  for (i = 0; i < m->size; i++)
    for (j = 0; j < m->size; j++)
      m_off.a[i][j] = m->a[i][j] + m->a[i][j] * RAND_VALUE(MAX_RAND_OFFSET);

  return m_off;
}

matr TransposeMatr( matr *m )
{
  matr tr = Unit(m->size);
  int i, j;

  for (i = 0; i < m->size; i++)
    for (j = 0; j < m->size; j++)
      tr.a[i][j] = m->a[j][i];

  return tr;
}

matr MakeSymPosMatr(int size)
{
  matr un = Unit(size), p;
  vec v, w = RandVec(size);
  matr vvt, vvt2;
  int i;

  for (i = 0; i < size; i++)
    while (w.v[i] == 0)
      w.v[i] = (int)RAND_VALUE(MAX_RAND_VALUE);

  v = VecMulNum(&w, 1 / NormEucl(&w));
  vvt = VecMulVec(&v, &v);
  vvt2 = MatrMulNum(&vvt, 2);

  p = MatrPlusMatr(&un, &vvt);

  FreeMatr(&un);
  FreeMatr(&vvt);
  FreeMatr(&vvt2);
  FreeVec(&v);
  FreeVec(&w);

  return p;
}

void CopyMatr(matr * m1, matr * m2)
{
  int i, j;

  for (i = 0; i < m1->size; i++)
    for (j = 0; j < m1->size; j++)
      m1->a[i][j] = m2->a[i][j];
}

double Det(matr * m)
{
  int i, j, k, pwr = 0;
  long double det = 1, treshold = 0.0001;
  matr m_t;

  if (m->size == 1)
    return m->a[0][0];
  else if (m->size == 2)
    return s_determ2x2(m);
  else if (s_triangle(m))
  {
    for (i = 0; i < m->size; i++)
      det *= m->a[i][i];

    return det;
  }
  else
  {
    m_t = Unit(m->size);
    CopyMatr(&m_t, m);

    for (j = 0; j < m_t.size; j++)
    {
      /* Swap rows to non-zero */
      if (m_t.a[j][j] == 0)
      {
        int old_pwr = pwr;

        for (i = j + 1; i < m_t.size; i++)
          if (m_t.a[i][j] != 0)
          {
            SwapPtr(&m_t.a[i], &m_t.a[j]);
            pwr++;
          }

        /* No swaps => no row with non-zero => det = 0 */
        if (pwr == old_pwr)
        {
          FreeMatr(&m_t);
          return 0;
        }
      }

      /* First Gauss step */
      for (i = j + 1; i < m_t.size; i++)
      {
        double prod = m_t.a[i][j] / m_t.a[j][j];

        for (k = j; k < m_t.size; k++)
          m_t.a[i][k] -= m_t.a[j][k] * prod;
      }

      if (s_triangle(&m_t))
        continue;
    }
  }

  for (i = 0; i < m_t.size; i++)
    det *= m_t.a[i][i];

  FreeMatr(&m_t);
  return det * pow(-1, pwr);
}

double GetCond(matr * m)
{
  matr m_t = InverseMatr(m);
  double cond = InfNorm(m) * InfNorm(&m_t);

  FreeMatr(&m_t);
  return cond;
}

double InfNorm(matr * m)
{
  double max = 0;
  int i, j;

  for (j = 0; j < m->size; j++)
    max += fabs(m->a[0][j]);

  for (i = 1; i < m->size; i++)
  {
    double sum = 0;

    for (j = 0; j < m->size; j++)
      sum += fabs(m->a[i][j]);

    if (sum > max)
      max = sum;
  }

  return max;
}

void PrintMatr(matr * m, char *name, FILE *f)
{
  int i, j;

  if (*name != 0)
    fprintf(f, "\n%s\n", name);

  for (i = 0; i < m->size; i++)
  {
    for (j = 0; j < m->size; j++)
    {
      if (m->a[i][j] >= 0)
        fprintf(f, " ");
      fprintf(f, "%7lf ", m->a[i][j]);
    }
    fprintf(f, "\n");
  }
}

void FreeMatr(matr * m)
{
  int i;

  for (i = 0; i < m->size; i++)
    free(m->a[i]);

  free(m->a);
}
