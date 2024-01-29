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
// spatial density function and its derivatives
// ******************************************************************

double fr(double r, double w[])
{
  double gamma = w[0], D = w[1], q = w[2], mag = w[3], sig = D * exp(gamma * mag);
  return (1 - pow(1 + r * r / sig, 1 - q)) / (2 * M_PI);
}

double dgamma_fr(double r, double w[])
{
  double gamma = w[0], D = w[1], q = w[2], mag = w[3], sig = D * exp(gamma * mag);
  return (1 -q) * pow(1 + r * r / sig, -q) * mag * r * r / sig / (2 * M_PI);
}

double dD_fr(double r, double w[])
{
  double gamma = w[0], D = w[1], q = w[2], mag = w[3], sig = D * exp(gamma * mag);
  return (1 - q) * pow(1 + r * r / sig, -q) / D * r * r / sig / (2 * M_PI);
}

double dq_fr(double r, double w[])
{
  double gamma = w[0], D = w[1], q = w[2], mag = w[3], sig = D * exp(gamma * mag);
  return pow(1 + r * r / sig, 1 - q) * log(1 + r * r / sig) / (2 * M_PI);
}

// ******************************************************************
#endif