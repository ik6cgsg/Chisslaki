#ifndef _INTEGR_H__
#define _INTEGR_H__
#pragma once

double Integrate( double (*f)( double ), double a, double b, int pNum );
double IntegrateRunge( double (*f)( double ), double a, double b, double eps, int *pNum );
double NewLeyb( double (*fInt)( double ), double a, double b );

#endif /* _INTEGR_H__ */
