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
// expected number of triggered events
// ******************************************************************

double kappafun(double m, double kparam[])
{
  double A = k_param[0], alpha = kparam[1];
  return A * exp(alpha * m);
}

NumericVector dkappafun(double m, double kparam[])
{
  double A = k_param[0], alpha = kparam[1];
  NumericVector out(3);
  out[0] = A * exp(alpha * m);
  // d A
  out[1] = out[0] / A;
  // d alpha
  out[2] = out[0] * m;
  return out;
}

// ******************************************************************
// temporal density function and its derivatives
// ******************************************************************

double gfun(double t, double gparam[])
{
  double c = gparam[0], p = gparam[1];
  if (t <= 0)
    return 0;
  else
    return (p - 1) / c * pow(1 + t / c, - p);
}

NumericVector dgfun(double t, double gparam[])
{
  double c = gparam[0], p = gparam[1];
  NumericVector out(3);
  out[0] = (p - 1) / c * pow(1 + t / c, - p);
  // d c
  out[1] = out[0] * (-1 / c - p / (c + t) + p / c);
  // d p
  out[2] = out[0] * (1 / (p - 1) - log(1 + t / c));
  return out;
}

double g1i(double t, double c, double p)
{
  return 1 - pow(1 + t / c, 1 - p);
}

NumericVector dgifun(double t, double c, double p)
{
  NumericVector out(3);
  out[0] = 1 - pow(1 + t / c, 1 - p);
  // d c
  out[1] = - (1 - out[0]) * (1 - p) * (1 / (c + t) - 1 / c);
  // d p
  out[2] = - (1 - out[0]) * (log(c) - log(c + t));
  return out;
}

// ******************************************************************
// spatial density function and its derivatives
// ******************************************************************

double f1(double r2, double sig, double q)
{
  return (q - 1) / (sig * M_PI) * pow(1 + r2 / sig, - q);
}

NumericVector dffun(double r2, double sig, double q)
{
  NumericVector out(3);
  out[0] = (q - 1) / (sig * M_PI) * pow(1 + r2 / sig, - q);
  // d sig
  out[1] = out[0] * (-1 + q * r2  / (r2 + sig)) / sig;
  // d q
  out[2] = out[0] * (1 / (q - 1) - log(1 + r2 / sig));
  return out;
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

NumericVector dfrifun(double r, double sig, double q)
{
  NumericVector out(3);
  double r2 = r * r / sig;
  double v = pow(1 + r2, 1 - q) / (2 * M_PI);
  out[0] = 1 / (2 * M_PI) - v;
  // d sig
  out[1] = (1 - q) * v / (1 + r2) * r2 / sig;
  // d q
  out[2] =  v * log(1 + r2);
  return out;
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
