#include "vec.h"

vec Zeros( int length )
{
  vec v;
  int i;

  v.v = malloc(sizeof(double) * length);
  v.length = length;

  for (i = 0; i < length; i++)
    v.v[i] = 0;

  return v;
}

vec Ones( int length )
{
  vec v;
  int i;

  v.v = malloc(sizeof(double) * length);
  v.length = length;

  for (i = 0; i < length; i++)
    v.v[i] = 1;

  return v;
}

vec RandVec( int length )
{
  vec v;
  int i;

  v.v = malloc(sizeof(double) * length);
  v.length = length;

  for (i = 0; i < length; i++)
    v.v[i] = (int)RAND_VALUE(MAX_RAND_VALUE);

  return v;
}

vec FullRandVec(int length)
{
  vec v = Zeros(length);
  int i;

  for (i = 0; i < length; i++)
    while (v.v[i] == 0)
      v.v[i] = (int)RAND_VALUE(MAX_RAND_VALUE);

  return v;
}

void CopyVec( vec * v1, vec * v2 )
{
  int i;

  for (i = 0; i < v1->length; i++)
    v1->v[i] = v2->v[i];
}

vec VecMulNum( vec * v, double n )
{
  int i;
  vec vn = Zeros(v->length);

  for (i = 0; i < v->length; i++)
    vn.v[i] = v->v[i] * n;

  return vn;
}

vec VecPlusVec( vec * v1, vec * v2 )
{
  int i;
  vec v = Zeros(v1->length);

  for (i = 0; i < v1->length; i++)
    v.v[i] = v1->v[i] + v2->v[i];

  return v;
}

vec VecMinusVec( vec * v1, vec * v2 )
{
  int i;
  vec v = Zeros(v1->length);

  for (i = 0; i < v1->length; i++)
    v.v[i] = v1->v[i] - v2->v[i];

  return v;
}

matr VecMulVec( vec * v1, vec * v2 )
{
  matr prod = Unit(v1->length);
  int i, j;

  for (i = 0; i < v1->length; i++)
    for (j = 0; j < v1->length; j++)
      prod.a[i][j] = v1->v[i] * v2->v[j];

  return prod;
}

double Scalar(vec * v1, vec * v2)
{
  int i;
  double sum = 0;

  for (i = 0; i < v1->length; i++)
    sum += v1->v[i] * v2->v[i];

  return sum;
}

vec MatrMulVec(matr * m, vec * v)
{
  vec prod = Zeros(m->size);
  int i;

  for (i = 0; i < m->size; i++)
  {
    double sum = 0;
    int k;

    for (k = 0; k < m->size; k++)
      sum += m->a[i][k] * v->v[k];

    prod.v[i] = sum;
  }

  return prod;
}

vec OffsetVec( vec * v )
{
  vec v_off = Zeros(v->length);
  int i;

  for (i = 0; i < v->length; i++)
    v_off.v[i] = v->v[i] + v->v[i] * RAND_VALUE(MAX_RAND_OFFSET);

  return v_off;
}

vec GetDiag( matr * m )
{
  vec diag = Zeros(m->size);
  int i;

  for (i = 0; i < m->size; i++)
    diag.v[i] = m->a[i][i];

  return diag;
}

matr Householder( int size, vec *w )
{
  matr un = Unit(size), p;
  matr vvt, vvt2;

  vvt = VecMulVec(w, w);
  vvt2 = MatrMulNum(&vvt, 2);

  p = MatrMinusMatr(&un, &vvt2);

  FreeMatr(&un);
  FreeMatr(&vvt);
  FreeMatr(&vvt2);

  return p;
}

double NormInf( vec * v )
{
  double max = fabs(v->v[0]);
  int i;

  for (i = 1; i < v->length; i++)
    if (fabs(v->v[i]) > max)
      max = fabs(v->v[i]);

  return max;
}

double NormEucl( vec *v )
{
  int i;
  double sum = 0;

  for (i = 0; i < v->length; i++)
    sum += pow(v->v[i], 2);

  return sqrt(sum);
}

void PrintVec( vec * v, char *name, FILE *f )
{
  int i;

  if (*name != 0)
    fprintf(f, "\n%s\n", name);

  for (i = 0; i < v->length; i++)
  {
    if (v->v[i] > 0)
      fprintf(f, " ");
    fprintf(f, "%.10lf\n", v->v[i]);
  }
}

void FreeVec( vec * v )
{
  free(v->v);
}
