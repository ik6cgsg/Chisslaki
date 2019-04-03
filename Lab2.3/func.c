#include <stdio.h>
#include <math.h>

#include "func.h"

double func1( double x )
{
  return exp(x);
} /* End of first function */

double func2( double x )
{
  return pow(exp(x), 5);
} /* End of second function */

double func1deriv4( double x )
{
  return exp(x);
} /* End of 'func1deriv4' function */

double func2deriv4( double x )
{
  return 625 * exp(5 * x);
} /* End of 'func2deriv4' function */

double func1int( double x )
{
  return exp(x);
} /* End of 'func1int' fucntion */

double func2int( double x )
{
  return exp(5 * x) / 5;
} /* End of 'func2int' function */
