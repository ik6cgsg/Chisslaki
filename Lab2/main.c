#include <stdio.h>

#include "decompos.h"

#define SHOW_LU(M, ML, MU) \
LUDecompose(&(M), &(ML), &(MU)); \
                                 \
PrintMatr(&(ML), "Left");        \
PrintMatr(&(MU), "Right");

#define SIZE 10

void main( void )
{
  matr m_hilb, m_rand, m_left, m_upper, m_hilb_off, m_rand_off, m_tmp, q, r;
  vec b, b1, x, x1, x_off, x1_off, b_off, b1_off, v_tmp, ax;
  double cond, a[9] = {12, -51, 4, 6, 167, -68, -4, 24, -41};
  FILE *f;
  vec (*solve)( matr*, vec* ) = SolveLU/*SolveQR*/;

  _CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );

#if 0
  m_hilb = RandMatr(3);
  PrintMatr(&m_hilb);

  b.v[0] = 1;
  b.v[1] = 2;
  b.v[2] = 3;

  printf("\n");
  PrintVec(&b);

  printf("\ndet = %g\n", Det(&m_hilb));

  m_left = Unit(3);
  m_upper = Unit(3);

  if (LUDecompose(&m_hilb, &m_left, &m_upper))
  {
    printf("\nLeft\n");
    PrintMatr(&m_left);
    
    printf("\nUpper\n");
    PrintMatr(&m_upper);
  }

  x = SolveLU(&m_hilb, &b);
  printf("\nX\n");
  PrintVec(&x);
#endif
#if 0
  //m_rand = RandMatr(10);
  m_rand = MakeMatr(3, a);
  b = RandVec(3);
  b.v[1] = 2;

  PrintMatr(&m_rand, "Rand", stdout);
  PrintVec(&b, "B", stdout);

  printf("\ndet = %g\n", Det(&m_rand));
  printf("\ncond(A) = %g\n", GetCond(&m_rand));

  x1 = solve(&m_rand, &b);

  //PrintVec(&x, "LU sol", stdout);
  PrintVec(&x1, "QR sol", stdout);

  /*
  if (QRDecompose(&m_rand, &q, &r))
  {
    PrintMatr(&q, "Q", stdout);
    PrintMatr(&r, "R", stdout);

    m_tmp = MatrMulMatr(&q, &r);
    PrintMatr(&m_tmp, "Rand again", stdout);
  }

  FreeMatr(&m_rand);
  FreeMatr(&q);
  FreeMatr(&r);
  FreeMatr(&m_tmp);*/

#endif
#if 1
  /* Hilbert */

  /// a
  printf("-------------------------------------------------------\nA)");

  m_hilb = Hilb(SIZE);
  x = RandVec(SIZE);

  //x = VecMulNum(&x, 27001);
  b = /*RandVec(SIZE)*/MatrMulVec(&m_hilb, &x);
  cond = GetCond(&m_hilb);

  PrintMatr(&m_hilb, "Hilbert", stdout);
  PrintVec(&x, "Solution", stdout);
  PrintVec(&b, "B = Hilbert * Solution", stdout);

  printf("\ndet = %g\n", Det(&m_hilb));
  printf("\ncond(A) = %g\n", cond);

  x = solve(&m_hilb, &b);

  PrintVec(&x, "Program solution", stdout);

  ax = MatrMulVec(&m_hilb, &x);
  v_tmp = VecMinusVec(&ax, &b);
  PrintVec(&v_tmp, "Discrepancy", stdout);

  FreeVec(&ax);
  FreeVec(&v_tmp);

  /// b
  printf("-------------------------------------------------------\nB)");

  m_hilb_off = OffsetMatr(&m_hilb);

  PrintMatr(&m_hilb_off, "Hilbert + offset", stdout);

  x_off = solve(&m_hilb_off, &b);

  PrintVec(&x_off, "Offset solution", stdout);

  b_off = VecMinusVec(&x_off, &x);

  m_tmp = MatrMinusMatr(&m_hilb_off, &m_hilb);

  printf("\n||dX||/||X|| = %.10lf <= %g * ||dA||/||A|| = %g * %lf\n", Norm(&b_off) / Norm(&x_off), cond, cond, InfNorm(&m_tmp) / InfNorm(&m_hilb));

  /// c
  printf("-------------------------------------------------------\nC)");

  FreeVec(&b_off);
  b_off = OffsetVec(&b);

  FreeVec(&x_off);

  x_off = solve(&m_hilb, &b_off);

  PrintVec(&b_off, "B + offset", stdout);
  PrintVec(&x_off, "Offset solution", stdout);

  v_tmp = VecMinusVec(&x_off, &x);

  FreeVec(&x_off);
  x_off = VecMinusVec(&b_off, &b);

  printf("\n||dX||/||X|| = %lf  ||dB||/||B|| = %lf\n", Norm(&v_tmp) / Norm(&x), Norm(&x_off) / Norm(&b));

  printf("\n\n\n");

  /* Rand */

  /// a
  printf("-------------------------------------------------------\nA)");

  m_rand = RandMatr(SIZE);
  b1 = RandVec(SIZE);
  cond = GetCond(&m_rand);

  PrintMatr(&m_rand, "Random", stdout);
  PrintVec(&b1, "B", stdout);

  printf("\ndet = %g\n", Det(&m_rand));
  printf("\ncond(A) = %g\n", cond);

  x1 = solve(&m_rand, &b1);

  PrintVec(&x1, "Solution", stdout);

  ax = MatrMulVec(&m_rand, &x1);
  v_tmp = VecMinusVec(&ax, &b1);
  PrintVec(&v_tmp, "Eps", stdout);

  FreeVec(&ax);
  FreeVec(&v_tmp);

  /// b
  printf("-------------------------------------------------------\nB)");

  m_rand_off = OffsetMatr(&m_rand);

  PrintMatr(&m_rand_off, "Random + offset", stdout);

  x1_off = solve(&m_rand_off, &b1);

  PrintVec(&x1_off, "Offset solution", stdout);

  b1_off = VecMinusVec(&x1_off, &x1);

  FreeMatr(&m_tmp);
  m_tmp = MatrMinusMatr(&m_rand_off, &m_rand);

  printf("\n||dX||/||X|| = %lf <= %g * ||dA||/||A|| = %g * %lf\n", Norm(&b1_off) / Norm(&x1), cond, cond, InfNorm(&m_tmp) / InfNorm(&m_rand));

  /// c
  printf("-------------------------------------------------------\nC)");

  FreeVec(&b1_off);
  b1_off = OffsetVec(&b1);

  FreeVec(&x1_off);

  x1_off = solve(&m_rand, &b1_off);

  PrintVec(&b1_off, "B + offset", stdout);
  PrintVec(&x1_off, "Offset solution", stdout);

  //FreeVec(&v_tmp);
  v_tmp = VecMinusVec(&x1_off, &x1);

  FreeVec(&x1_off);
  x1_off = VecMinusVec(&b1_off, &b1);

  printf("\n||dX||/||X|| = %lf  ||dB||/||B|| = %lf\n", Norm(&v_tmp) / Norm(&x1), Norm(&x1_off) / Norm(&b1));

  /// print matrixes and rights parts

  f = fopen("randmatr", "w");
  PrintMatr(&m_rand, "", f);
  fclose(f);

  f = fopen("randmatroff", "w");
  PrintMatr(&m_rand_off, "", f);
  fclose(f);

  f = fopen("hilboff", "w");
  PrintMatr(&m_hilb_off, "", f);
  fclose(f);

  f = fopen("bhilb", "w");
  PrintVec(&b, "", f);
  fclose(f);

  f = fopen("bhilboff", "w");
  PrintVec(&b_off, "", f);
  fclose(f);

  f = fopen("brand", "w");
  PrintVec(&b1, "", f);
  fclose(f);

  f = fopen("brandoff", "w");
  PrintVec(&b1_off, "", f);
  fclose(f);

  /// free

  FreeVec(&b);
  FreeVec(&b1);
  FreeVec(&x);
  FreeVec(&x1);
  FreeVec(&x_off);
  FreeVec(&x1_off);
  FreeVec(&b_off);
  FreeVec(&b1_off);
  FreeVec(&v_tmp);

  FreeMatr(&m_hilb);
  FreeMatr(&m_hilb_off);
  FreeMatr(&m_rand);
  FreeMatr(&m_rand_off);
  //FreeMatr(&m_left);
  //FreeMatr(&m_upper);
  FreeMatr(&m_tmp);

#endif
}
