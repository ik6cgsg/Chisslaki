#include <stdio.h>

#include "interpol.h"

#define GRID_SIZE 10

double func1( double x )
{
  return exp(x);
}

double func2( double x )
{
  return exp(x) / (1 - x * x);
}

double maxDelta( double (*f)(double), grid *gr )
{
  uint i;
  double max = 0;

  for (i = 1; i < gr->size; i++)
  {
    double sub;
    double x = (gr->nodes[i - 1] + gr->nodes[i]) / 2;

    sub = fabs(InterPolLagrange(f, gr, x) - f(x));

    if (sub > max)
      max = sub;
  }

  return max;
}

void generateFile( char const *name, int points, int nodes, grid_type gt )
{
  FILE *fout;
  grid gr, gr_y;
  int i;

  if ((fout = fopen(name, "w")) == NULL)
    return 1;

  GridInit(&gr, points);
  GridInit(&gr_y, gr.size);

  GridUniform(&gr, -1, 1);
  GenerateY(&gr, &gr_y, func1, nodes, gt);

  fprintf(fout, "%li\n", points);

  for (i = 0; i < gr.size; i++)
  {
    fprintf(fout, "%lf ", gr.nodes[i]);
  }
  fprintf(fout, "\n");

  for (i = 0; i < gr.size; i++)
  {
    fprintf(fout, "%lf ", gr_y.nodes[i]);
  }
  fprintf(fout, "\n");

  GridUniform(&gr, -1 + 0.08, 1 - 0.08);
  for (i = 0; i < gr.size; i++)
  {
    fprintf(fout, "%lf ", gr.nodes[i]);
  }
  fprintf(fout, "\n");

  GenerateY(&gr, &gr_y, func2, nodes, gt);
  for (i = 0; i < gr.size; i++)
  {
    fprintf(fout, "%lf ", gr_y.nodes[i]);
  }

  GridFree(&gr);
  GridFree(&gr_y);

  fclose(fout);
}

int main( void )
{
  grid gr, gr_uni, gr_my;
  uint i, pnts = 2000;
  double (*f)(double) = func1;
  double d = 0.2;

  /* Generate file with grids */
#if 1
  generateFile("f.out", pnts, 3, GRID_UNIFORM);
  generateFile("s.out", pnts, 5, GRID_UNIFORM);
  generateFile("t.out", pnts, 9, GRID_UNIFORM);

  generateFile("ff.out", pnts, 3, GRID_CHEBYSHEV);
  generateFile("ss.out", pnts, 5, GRID_CHEBYSHEV);
  generateFile("tt.out", pnts, 9, GRID_CHEBYSHEV);

  generateFile("fff.out", pnts, 3, GRID_USER);
  generateFile("sss.out", pnts, 5, GRID_USER);
  generateFile("ttt.out", pnts, 9, GRID_USER);

  getchar();
#endif
  /* Experiments */
#if 0
  GridInit(&gr, GRID_SIZE);

  printf("exp(x)\n\n1) grid size: 10\n");
  GridUniform(&gr, -1, 1);
  printf("Uniform grid delta: %g\n", maxDelta(f, &gr));
  GridChebyshev(&gr, -1, 1);
  printf("Chebyshev grid delta: %g\n", maxDelta(f, &gr));
  GridHilb(&gr, -1, 1);
  printf("Hilb grid delta: %g\n\n", maxDelta(f, &gr));

  printf("2)\nnodes number | delta (Chebyshev) | delta (Uni) | delta (my)\n");
  printf("-----------------------------------------------------------\n");
  for (i = 2; i < 21; i++)
  {
    GridFree(&gr);
    GridInit(&gr, i);
    GridInit(&gr_uni, i);
    GridInit(&gr_my, i);

    GridChebyshev(&gr, -1, 1);
    GridUniform(&gr_uni, -1, 1);
    GridHilb(&gr_my, -1, 1);

    printf("%2i           | %11g       | %11g | %11g\n", i, maxDelta(f, &gr), maxDelta(f, &gr_uni), maxDelta(f, &gr_my));

    GridFree(&gr_uni);
    GridFree(&gr_my);
  }

  GridFree(&gr);

  printf("\n\n");

  f = func2;

  GridInit(&gr, GRID_SIZE);

  printf("exp(x) / (1 - x * x)\n\n1) grid size: 10 interval delta: 0.2\n");
  GridUniform(&gr, -1 + d, 1 - d);
  printf("Uniform grid delta: %g\n", maxDelta(f, &gr));
  GridChebyshev(&gr, -1 + d, 1 - d);
  printf("Chebyshev grid delta: %g\n", maxDelta(f, &gr));
  GridHilb(&gr, -1 + d, 1 - d);
  printf("Hilb grid delta: %g\n\n", maxDelta(f, &gr));

  printf("2)\nnodes number | delta (Chebyshev) | delta (Uni) | delta (my)\n");
  printf("-----------------------------------------------------------\n");
  for (i = 2; i < 21; i++)
  {
    GridFree(&gr);
    GridInit(&gr, i);
    GridInit(&gr_uni, i);
    GridInit(&gr_my, i);

    GridChebyshev(&gr, -1 + d, 1 - d);
    GridUniform(&gr_uni, -1 + d, 1 - d);
    GridHilb(&gr_my, -1 + d, 1 - d);
    printf("%2i           | %11g       | %11g | %11g\n", i, maxDelta(f, &gr), maxDelta(f, &gr_uni), maxDelta(f, &gr_my));

    GridFree(&gr_uni);
    GridFree(&gr_my);
  }

  GridFree(&gr);

  printf("\n3)\ninterval delta | delta (Chebyshev) | delta (Uni) | delta (my) (grid size: 10)\n");
  printf("-----------------------------------------------------------\n");

  GridInit(&gr, 10);
  GridInit(&gr_uni, 10);
  GridInit(&gr_my, 10);

  for (i = 1; i < 10; i++)
  {
    d = pow(0.5, 2 * i);

    GridChebyshev(&gr, -1 + d, 1 - d);
    GridUniform(&gr_uni, -1 + d, 1 - d);
    GridHilb(&gr_my, -1 + d, 1 - d);
    printf("%.8lf     | %11g       | %9g   | %g\n", d, maxDelta(f, &gr), maxDelta(f, &gr_uni), maxDelta(f, &gr_my));

  }

  GridFree(&gr_uni);
  GridFree(&gr_my);
  GridFree(&gr);
#endif
  /*
  printf("Uniform grid:\n");
  GridUniform(&gr, -1, 1);
  for (i = 0; i < GRID_SIZE; i++)
    printf("%g ", gr.nodes[i]);
  printf("\n");

  printf("Lagrange values:\n");
  for (i = 0; i < GRID_SIZE; i++)
    printf("%g ", InterPolLagrange(f, &gr, gr.nodes[i]));
  printf("\n");

  printf("Original values:\n");
  for (i = 0; i < GRID_SIZE; i++)
    printf("%g ", f(gr.nodes[i]));
  printf("\n\n");

  printf("Chebyshev grid:\n");
  GridChebyshev(&gr, -1, 1);
  for (i = 0; i < GRID_SIZE; i++)
    printf("%g ", gr.nodes[i]);
  printf("\n");

  printf("Lagrange values:\n");
  for (i = 0; i < GRID_SIZE; i++)
    printf("%g ", InterPolLagrange(f, &gr, gr.nodes[i]));
  printf("\n");

  printf("Original values:\n");
  for (i = 0; i < GRID_SIZE; i++)
    printf("%g ", f(gr.nodes[i]));
  printf("\n\n");

  printf("Hilb grid:\n");
  GridHilb(&gr, -1, 1);
  for (i = 0; i < GRID_SIZE; i++)
    printf("%g ", gr.nodes[i]);
  printf("\n");

  printf("Lagrange values:\n");
  for (i = 0; i < GRID_SIZE; i++)
    printf("%g ", InterPolLagrange(f, &gr, gr.nodes[i]));
  printf("\n");

  printf("Original values:\n");
  for (i = 0; i < GRID_SIZE; i++)
    printf("%g ", f(gr.nodes[i]));
  printf("\n\n");
  */

  getchar();

  return 0;
}
