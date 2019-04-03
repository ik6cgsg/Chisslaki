#include <stdio.h>

#include "func.h"

double func1( double x )
{
  return exp(x);
} /* End of first function */

double func2( double x )
{
  return exp(x) / (1 - x * x);
} /* End of second function */

double secDerivFunc2( double x )
{
  return exp(x) * (pow(x, 4) - 4 * pow(x, 3) + 4 * pow(x, 2) + 4 * x + 3) / pow(1 - x * x, 3);
} /* End of 'secDerivFunc2' function */

double maxDelta( double (*f)(double), spline *sp )
{
  uint i;
  double max = 0;

  for (i = 1; i < sp->gr.size; i++)
  {
    double sub = 0;
    double x = (sp->gr.nodes[i - 1] + sp->gr.nodes[i]) / 2;

    sub = fabs(GetValue(sp, x) - f(x));

    if (sub > max)
      max = sub;
  }

  return max;
} /* End of 'maxDelta' function */

void GenerateY( grid const *x, grid *y_nat, grid *y_dev, double (*f)(double), double dev1, double dev2, uint nodes_num )
{
  uint i;
  spline sp;

  SplineInit(&sp, nodes_num, x->nodes[0], x->nodes[x->size - 1]);
  BuildNaturalSpline(&sp, f);

  for (i = 0; i < x->size; i++)
    y_nat->nodes[i] = GetValue(&sp, x->nodes[i]);

  BuildSecDerivSpline(&sp, f, dev1, dev2);

  for (i = 0; i < x->size; i++)
    y_dev->nodes[i] = GetValue(&sp, x->nodes[i]);

  SplineFree(&sp);
} /* End of 'GenerateY' function */

void generateFile( char const *name, int points, int nodes )
{
  FILE *fout;
  grid gr, gr_y, gr_y2;
  uint i;

  if ((fout = fopen(name, "w")) == NULL)
    return;

  GridInit(&gr, points);
  GridInit(&gr_y, gr.size);
  GridInit(&gr_y2, gr.size);

  GridUniform(&gr, -1, 1);
  GenerateY(&gr, &gr_y, &gr_y2, func1, func1(gr.nodes[0]), func1(gr.nodes[gr.size - 1]), nodes);

  fprintf(fout, "%li\n", points);

  for (i = 0; i < gr.size; i++)
    fprintf(fout, "%lf ", gr.nodes[i]);
  fprintf(fout, "\n");

  for (i = 0; i < gr.size; i++)
    fprintf(fout, "%lf ", gr_y.nodes[i]);
  fprintf(fout, "\n");

  for (i = 0; i < gr.size; i++)
    fprintf(fout, "%lf ", gr_y2.nodes[i]);
  fprintf(fout, "\n");

  GridUniform(&gr, -1 + 0.05, 1 - 0.05);

  for (i = 0; i < gr.size; i++)
    fprintf(fout, "%lf ", gr.nodes[i]);
  fprintf(fout, "\n");

  GenerateY(&gr, &gr_y, &gr_y2, func2, secDerivFunc2(gr.nodes[0]), secDerivFunc2(gr.nodes[gr.size - 1]), nodes);

  for (i = 0; i < gr.size; i++)
    fprintf(fout, "%lf ", gr_y.nodes[i]);
  fprintf(fout, "\n");

  for (i = 0; i < gr.size; i++)
    fprintf(fout, "%lf ", gr_y2.nodes[i]);
  fprintf(fout, "\n");

  GridFree(&gr);
  GridFree(&gr_y);

  fclose(fout);
} /* End of 'generateFile' function */
