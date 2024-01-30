#ifndef FUNCS_H
#define FUNCS_H
// ******************************************************************

#include <Rcpp.h>
using namespace Rcpp;

// ******************************************************************
// distance and norm functions
// ******************************************************************

inline int sgn(double x)
{
  return (x < 0) ? -1 : 1;
}

inline
double dist(double x1, double y1, double x2, double y2)
{
  return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

inline
double dist2(double x1, double y1, double x2, double y2)
{
  return ((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

double norm(double *x, int dim)
{
  double sum = 0;
  for (int i = 0; i < dim; i++)
    sum += x[i] * x[i];
  return sqrt(sum);
}

// ******************************************************************
// temporal density function and its derivatives
// ******************************************************************

double g1(double t, double c, double p)
{
  if (t <= 0)
    return 0;
  else
    return (p - 1) / c * pow(1 + t / c, - p);
}

double g1i(double a, double b, double c, double p)
{
  if (a == 0)
    return 1 - pow(1 + b / c, 1 - p);
  else
    return pow(1 + a / c, 1 - p) - pow(1 + b / c, 1 - p);
}

// ******************************************************************
// spatial density function and its derivatives
// ******************************************************************

double f1(double r2, double sig, double q)
{
  return (q - 1) / (sig * M_PI) * pow(1 + r2 / sig, - q);
}

double dq_f1(double r2, double sig, double q)
{
  return (1/(q - 1) - log(1 + r2/sig));
}

double dsig_f1(double r2, double sig, double q)
{
  return  (-1 + q * r2 /(r2 + sig)) / sig;
}

double f1r(double r, double w[])
{
  double sig = w[0], q = w[1];
  return (1 - pow(1 + r * r / sig, 1 - q)) / (2 * M_PI);
}

double dq_f1r(double r, double w[])
{
  double sig = w[0], q = w[1];
  return pow(1 + r * r / sig, 1 - q) * log(1 + r * r / sig) / (2 * M_PI);
}

double dsig_f1r(double r, double w[])
{
  double sig = w[0], q = w[1];
  return (1 - q) * pow(1 + r * r / sig, -q) * r * r / (sig * sig) / (2 * M_PI);
}

// ******************************************************************
#endif
