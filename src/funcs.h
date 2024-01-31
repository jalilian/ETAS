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

double* dgfun(double t, double c, double p)
{
  static double out[3];
  out[0] = g1(t, c, p);
  // d c
  out[1] = out[0] * (-1 / c - p / (c + t) + p / c);
  // d p
  out[2] = out[0] * (1 / (p - 1) - log(1 + t / c));
  return out;
}

// d g1 / g1
double dc_g1(double t, double c, double p)
{
  return -1 / c - p / (c + t) + p / c;
}

double dp_g1(double t, double c, double p)
{
  return 1 / (p - 1) - log(1 + t / c);
}

double g1i(double t, double c, double p)
{
  return 1 - pow(1 + t / c, 1 - p);
}

double dc_g1i(double t, double c, double p)
{
  return (1 - p) * t / (c * c) * pow(1 + t / c, -p);
}

double dp_g1i(double t, double c, double p)
{
  return pow(1 + t / c, 1 - p) * log(1 + t / c);
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
