#include "decompos.h"

int LUDecompose(matr * orig, matr * l, matr * u)
{
  int i, j, k;
  matr org = Unit(orig->size);
  CopyMatr(&org, orig);

  for (i = 0; i < org.size; i++)
    for (j = 0; j < org.size; j++)
      if (i < j)
      {
        double sum = 0;

        for (k = 0; k <= i - 1; k++)
          sum += l->a[i][k] * u->a[k][j];

        u->a[i][j] = org.a[i][j] - sum;
      }
      else if (i == j)
      {
        double sum = 0;
        int e = org.size - 1;

        for (k = 0; k <= i - 1; k++)
          sum += l->a[i][k] * u->a[k][i];

        u->a[i][i] = org.a[i][i] - sum;

        /* swap */
        if (u->a[i][i] == 0)
        {
          while (u->a[i][i] == 0 && e > i)
          {
            SwapPtr(&org.a[i], &org.a[e]);
            u->a[i][i] = org.a[i][i];
            e--;
          }

          if (e == i)  // det == 0
          {
            FreeMatr(&org);
            return 0;
          }
          else
            j = 0;  // restart
        }
      }
      else
      {
        double sum = 0;

        for (k = 0; k <= j - 1; k++)
          sum += l->a[i][k] * u->a[k][j];

        l->a[i][j] = (org.a[i][j] - sum) / u->a[j][j];
      }

  FreeMatr(&org);
  return 1;
}

static void s_putOnes(matr *m, int index)
{
  int i, j;

  for (i = 0; i <= index; i++)
  {
    int k;

    for (j = i; j < m->size; j++)
      m->a[i][j] = j == i ? 1 : 0;

    for (k = i + 1; k < m->size; k++)
      m->a[k][i] = 0;
  }
}

static vec s_getCol(matr *m, int index)
{
  int i, j;
  vec col = Zeros(m->size);

  for (i = 0; i < m->size; i++)
    col.v[i] = m->a[i][index];

  return col;
}

int QRDecompose(matr * orig, matr * q, matr * r)
{
  int i;
  matr q_t = Unit(orig->size);
  matr a = Unit(orig->size);

  if (Det(orig) == 0)
    return 0;

  CopyMatr(&a, orig);

  for (i = 0; i < orig->size - 1; i++)
  {
    vec e_i = Zeros(orig->size);
    vec x;
    vec u, v;
    matr tmp, q_i, unit = Unit(orig->size);
    double norm;

    s_putOnes(&a, i - 1);

    x = s_getCol(&a, i);
    norm = Len(&x);

    e_i.v[i] = x.v[i] > 0 ? -1 * norm : norm;

    u = VecMinusVec(&x, &e_i);
    norm = Len(&u);
    v = VecMulNum(&u, 1 / norm);
    q_i = VecMulVec(&v, &v);
    tmp = MatrMulNum(&q_i, 2);
    FreeMatr(&q_i);
    q_i = MatrMinusMatr(&unit, &tmp);
    FreeMatr(&tmp);
    tmp = MatrMulMatr(&q_i, &a, i);
    CopyMatr(&a, &tmp);

    FreeMatr(&unit);
    FreeMatr(&tmp);

    tmp = MatrMulMatr(&q_i, &q_t, 0);
    CopyMatr(&q_t, &tmp);

    FreeMatr(&tmp);
    FreeMatr(&q_i);
    FreeVec(&e_i);
    FreeVec(&u);
    FreeVec(&v);
    FreeVec(&x);
  }

  *q = TransposeMatr(&q_t);
  *r = MatrMulMatr(&q_t, orig, 0);

  FreeMatr(&q_t);
  FreeMatr(&a);

  return 1;
}

vec SolveLU(matr * m, vec * b)
{
  vec x = Zeros(m->size), y = Zeros(m->size);
  int i, j;
  matr m_l = Unit(m->size), m_u = Unit(m->size);

  LUDecompose(m, &m_l, &m_u);

  /* L(UX) = B => LY = B */
  for (i = 0; i < m->size; i++)
  {
    double sum = 0;

    for (j = 0; j < i; j++)
      sum += m_l.a[i][j] * y.v[j];

    y.v[i] = (b->v[i] - sum);
  }

  /* UX = Y */
  for (i = m->size - 1; i >= 0; i--)
  {
    double sum = 0;

    for (j = i + 1; j < m->size; j++)
      sum += m_u.a[i][j] * x.v[j];

    x.v[i] = (y.v[i] - sum) / m_u.a[i][i];
  }

  FreeMatr(&m_l);
  FreeMatr(&m_u);
  FreeVec(&y);

  return x;
}

vec SolveQR(matr * m, vec * b)
{
  vec x = Zeros(m->size), z;
  int i, j;
  matr q, r, tmp;

  if (!QRDecompose(m, &q, &r))
    return x;

  /*Z = Q^T x B*/
  tmp = TransposeMatr(&q);
  z = MatrMulVec(&tmp, b);

  /* RX = Z */
  for (i = m->size - 1; i >= 0; i--)
  {
    double sum = 0;

    for (j = i; j < m->size; j++)
      sum += r.a[i][j] * x.v[j];

    x.v[i] = (z.v[i] - sum) / r.a[i][i];
  }

  FreeMatr(&q);
  FreeMatr(&r);
  FreeMatr(&tmp);
  FreeVec(&z);

  return x;
}
