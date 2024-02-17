#ifndef FUNCS_H
#define FUNCS_H
// ******************************************************************

#include <array>
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
  double A = kparam[0], alpha = kparam[1];
  return A * exp(alpha * m);
}

std::array<double, 3> dkappafun(double m, double kparam[])
{
  double A = kparam[0], alpha = kparam[1];
  std::array<double, 3> out;
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

std::array<double, 3> dgfun(double t, double gparam[])
{
  double c = gparam[0], p = gparam[1];
  std::array<double, 3> out;
  out[0] = (p - 1) / c * pow(1 + t / c, - p);
  // d c
  out[1] = out[0] * (-1 / c - p / (c + t) + p / c);
  // d p
  out[2] = out[0] * (1 / (p - 1) - log(1 + t / c));
  return out;
}


double gfunint(double t, double gparam[])
{
  double c = gparam[0], p = gparam[1];
  return 1 - pow(1 + t / c, 1 - p);
}

std::array<double, 3> dgfunint(double t, double gparam[])
{
  double c = gparam[0], p = gparam[1];
  std::array<double, 3> out;
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

double ffun1(double r2, double m, double fparam[])
{
  double D = fparam[0], gamma = fparam[1], q = fparam[2];
  double sig = D * exp(gamma * m);
  return (q - 1) / (sig * M_PI) * pow(1 + r2 / sig, - q);
}

std::array<double, 4> dffun1(double r2, double m, double fparam[])
{
  double D = fparam[0], gamma = fparam[1], q = fparam[2];
  double sig = D * exp(gamma * m);
  std::array<double, 4> out;
  out[0] = (q - 1) / (sig * M_PI) * pow(1 + r2 / sig, - q);
  // d D
  out[1] = out[0] * (-1 + q * r2  / (r2 + sig)) / D;
  // d q
  out[2] = out[0] * (1 / (q - 1) - log(1 + r2 / sig));
  // d gamma
  out[3] = out[0] * (-1 + q * r2  / (r2 + sig)) * m;
  return out;
}

double ffunrint1(double r, double m, double fparam[])
{
  double D = fparam[0], gamma = fparam[1], q = fparam[2];
  double sig = D * exp(gamma * m);
  return (1 - pow(1 + r * r / sig, 1 - q)) / (2 * M_PI);
}

std::array<double, 4> dffunrint1(double r, double m, double fparam[])
{
  double D = fparam[0], gamma = fparam[1], q = fparam[2];
  double sig = D * exp(gamma * m);
  double r2 = r * r / sig;
  double v = pow(1 + r2, 1 - q) / (2 * M_PI);
  std::array<double, 4> out;
  out[0] = 1 / (2 * M_PI) - v;
  // d D
  out[1] = (1 - q) * v / (1 + r2) * r2 / D;
  // d q
  out[2] =  v * log(1 + r2);
  // d gamma
  out[3] = (1 - q) * v / (1 + r2) * r2 * m;
  return out;
}

// Gaussian density

double ffun2(double r2, double m, double fparam[])
{
  double D = fparam[0], gamma = fparam[1];
  double sig = D * exp(gamma * m);
  return exp(-r2 / (2 * sig * sig)) / (2 * M_PI * sig * sig);
}

std::array<double, 3> dffun2(double r2, double m, double fparam[])
{
  double D = fparam[0], gamma = fparam[1];
  double sig = D * exp(gamma * m);
  std::array<double, 3> out;
  out[0] = exp(-r2 / (2 * sig * sig)) / (2 * M_PI * sig * sig);
  // d D
  out[1] = out[0] * (r2 / (sig * sig) - 2) / D;
  // d gamma
  out[2] = out[0] * (r2 / (sig * sig) - 2) * m;
  return out;
}

double ffunrint2(double r, double m, double fparam[])
{
  double D = fparam[0], gamma = fparam[1];
  double sig = D * exp(gamma * m);
  return (1 - exp(-(r * r) / (2 * sig * sig))) / (2 * M_PI);
}

std::array<double, 3> dffunrint2(double r, double m, double fparam[])
{
  double D = fparam[0], gamma = fparam[1];
  double sig = D * exp(gamma * m);
  double r2 = r * r / (sig * sig);
  double v = exp(-r2 / 2) / (2 * M_PI);
  std::array<double, 3> out;
  out[0] = 1 / (2 * M_PI) - v;
  // d D
  out[1] = -v * r2 / D;
  // d gamma
  out[2] = -v * r2 * m;
  return out;
}

// ******************************************************************
#endif
