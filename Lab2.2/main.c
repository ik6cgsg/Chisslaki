#include <stdio.h>

#include "func.h"

void main( void )
{
  spline sp;
  uint i;
  double d = 0.1;

#if 1
  generateFile("a.out", 2000, 3);
  generateFile("b.out", 2000, 7);
#endif

#if 0
  printf("natural spline\n\nexp(x)\n\n");
  printf("nodes | delta\n------------------------\n");

  for (i = 3; i <= 20; i++)
  {
    SplineInit(&sp, i, -1, 1);

    BuildNaturalSpline(&sp, func1);
    printf("%2i    | %g\n", i, maxDelta(func1, &sp));

    SplineFree(&sp);
  }

  printf("\n\nexp(x) / (1 - x * x), delta = %g\n\n", d);
  printf("nodes | delta\n------------------------\n");

  for (i = 3; i <= 20; i++)
  {
    SplineInit(&sp, i, -1 + d, 1 - d);

    BuildNaturalSpline(&sp, func2);
    printf("%2i    | %g\n", i, maxDelta(func2, &sp));

    SplineFree(&sp);
  }

  printf("\n\n\nsecond derivative natural values spline\n\nexp(x)\n\n");
  printf("nodes | delta\n------------------------\n");

  for (i = 3; i <= 20; i++)
  {
    SplineInit(&sp, i, -1, 1);

    BuildSecDerivSpline(&sp, func1, func1(sp.gr.nodes[0]), func1(sp.gr.nodes[sp.gr.size - 1]));
    printf("%2i    | %g\n", i, maxDelta(func1, &sp));

    SplineFree(&sp);
  }

  printf("\n\nexp(x) / (1 - x * x), delta = %g\n\n", d);
  printf("nodes | delta\n------------------------\n");

  for (i = 3; i <= 20; i++)
  {
    SplineInit(&sp, i, -1 + d, 1 - d);

    BuildSecDerivSpline(&sp, func2, secDerivFunc2(sp.gr.nodes[0]), secDerivFunc2(sp.gr.nodes[sp.gr.size - 1]));
    printf("%2i    | %g\n", i, maxDelta(func2, &sp));

    SplineFree(&sp);
  }

  getch();

#endif

}
