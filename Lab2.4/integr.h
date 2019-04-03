#ifndef _INTEGR_H__
#define _INTEGR_H__
#pragma once

double Integrate38( double (*f)( double ), double a, double b, int pNum );
double IntegrateGauss( double (*f)( double ), double a, double b, int pNum );
double Integrate38ByEps( double (*f)( double ), double a, double b, double eps, int *pNum );
double IntegrateGaussByEps( double (*f)( double ), double a, double b, double eps, double res, int *pNum );
double NewLeib( double (*fInt)( double ), double a, double b );

#endif /* _INTEGR_H__ */
