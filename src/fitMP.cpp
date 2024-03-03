#include <climits>
#include <cmath>
#include <Rcpp.h>

#include "funcs.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

class modelhandler{
  private:
    int mver;
    double mu;
    double kparam[2];
    double gparam[2];
    double fparam[3];
  public:
    void set(int rmver, NumericVector param);
    double mufun(void);
    double kappafun0(double m);
    double gfun0(double t);
    double gfunint0(double t);
    double ffun0(double r2, double m);
    double ffunrint0(double r, double m);
};
void modelhandler::set(int rmver, NumericVector param)
{
  mver = rmver;

  mu = param[0];

  kparam[0] = param[1]; // A
  kparam[1] = param[3]; // alpha

  gparam[0] = param[2]; // c
  gparam[1] = param[4]; // p

  switch (mver)
  {
    case 1:
      fparam[0] = param[5]; // D
      fparam[1] = param[7]; // gamma
      fparam[2] = param[6]; // q
      break;
    case 2:
      fparam[0] = param[5]; // D
      fparam[1] = param[6]; // gamma
      break;
  }
}
double modelhandler::mufun(void)
{
  return mu;
}
double modelhandler::kappafun0(double m)
{
  return kappafun(m, kparam);
}
double modelhandler::gfun0(double t)
{
  return gfun(t, gparam);
}
double modelhandler::gfunint0(double t)
{
  return gfunint(t, gparam);
}
double modelhandler::ffun0(double r2, double m)
{
  double f = 0;
  switch (mver)
  {
    case 1:
      f = ffun1(r2, m, fparam);
      break;
    case 2:
      f = ffun2(r2, m, fparam);
      break;
  }
  return f;
}
double modelhandler::ffunrint0(double r, double m)
{
  double f = 0;
  switch (mver)
  {
    case 1:
      f = ffunrint1(r, m, fparam);
      break;
    case 2:
      f = ffunrint2(r, m, fparam);
      break;
  }
  return f;
}

// ******************************************************************
// the etas class
// ******************************************************************

class etas{
private:
  int N;
  NumericVector t;
  NumericVector x;
  NumericVector y;
  NumericVector m;
  NumericVector flag;
  NumericVector bk;
  NumericVector pb;
  NumericVector lam;
  int np;
  NumericVector px;
  NumericVector py;
  double tstart2;
  double tlength;
  double integ0;
  int ndiv;
  int mver;

public:
  void set(NumericMatrix revents,
           NumericMatrix rpoly,
           NumericVector tperiod,
           double rinteg0,
           int rndiv,
           int rmver);
  void paramhandler(NumericVector theta,
                    double *mu,
                    double *kparam,
                    double *gparam,
                    double *fparam);
  double mloglikj1(int j,
                   double mu,
                   double kparam[],
                   double gparam[],
                   double fparam[]);
  void mloglikj1Gr(int j,
                   double mu,
                   double kparam[],
                   double gparam[],
                   double fparam[],
                   double *fv,
                   double *df);
  double mloglikj2(int j,
                   double mu,
                   double kparam[],
                   double gparam[],
                   double fparam[]);
  void mloglikj2Gr(int j,
                   double mu,
                   double kparam[],
                   double gparam[],
                   double fparam[],
                   double *fv,
                   double *df);
  double mloglikj(int j,
                  double mu,
                  double kparam[],
                  double gparam[],
                  double fparam[]);
  void mloglikjGr(int j,
                  double mu,
                  double kparam[],
                  double gparam[],
                  double fparam[],
                  double *fv,
                  double *df);
  double mloglik(NumericVector theta);
  NumericVector mloglikGr(NumericVector theta);
  void linesearch(NumericVector xOld,
                  NumericVector h,
                  double *fv,
                  double *ram);
  List fitfun(NumericVector tht,
              NumericMatrix ihess,
              double eps,
              bool verbose);
  double mloglikMP(NumericVector theta,
                   int nthreads);
  NumericVector mloglikGrMP(NumericVector theta, int nthreads);
  void linesearchMP(NumericVector xOld,
                    NumericVector h,
                    double *fv,
                    double *ram,
                    int nthreads);
  List fitfunMP(NumericVector tht,
                NumericMatrix ihess,
                double eps,
                bool verbose,
                int nthreads);
  
};

// ******************************************************************
// set function
// ******************************************************************

void etas::set(NumericMatrix revents,
               NumericMatrix rpoly,
               NumericVector tperiod,
               double rinteg0,
               int rndiv,
               int rmver)
{
  N = revents.nrow();
  t = revents( _, 0);
  x = revents( _, 1);
  y = revents( _, 2);
  m = revents( _, 3);
  flag = revents( _, 4);
  bk = revents( _, 5);
  pb = revents( _, 6);
  lam = revents( _, 7);
  
  np = rpoly.nrow();
  px = rpoly( _, 0);
  py = rpoly( _, 1);
  
  tstart2 = tperiod[0];
  tlength = tperiod[1];
  
  integ0 = rinteg0;
  ndiv = rndiv;

  mver = rmver;
}

// ******************************************************************
// parameters of the model
// ******************************************************************

void etas::paramhandler(NumericVector theta,
                        double *mu,
                        double *kparam,
                        double *gparam,
                        double *fparam)
{
  *mu = theta[0] * theta[0];

  kparam[0] = theta[1] * theta[1]; // A
  kparam[1] = theta[3] * theta[3]; // alpha

  gparam[0] = theta[2] * theta[2]; // c
  gparam[1] = theta[4] * theta[4]; // p

  switch (mver)
  {
    case 1:
      fparam[0] = theta[5] * theta[5]; // D
      fparam[1] = theta[7] * theta[7]; // gamma
      fparam[2] = theta[6] * theta[6]; // q
      break;
    case 2:
      fparam[0] = theta[5] * theta[5]; // D
      fparam[1] = theta[6] * theta[6]; // gamma
      break;
  }
}

// ******************************************************************
// minus log likelihood function
// ******************************************************************

double etas::mloglikj1(int j,
                       double mu,
                       double kparam[],
                       double gparam[],
                       double fparam[])
{
  double sumpart = 0;
  if (flag[j] == 1)
  {
    double sumj = mu * bk[j];
    for (int i = 0; i < j; i++)
    {
      sumj += kappafun(m[i], kparam) *
        gfun(t[j] - t[i], gparam) *
        ffun1(dist2(x[j], y[j], x[i], y[i]), m[i], fparam);
    }

    sumpart = (sumj > 1.0e-25) ? log(sumj) : -100.0;
  }

  double gi = gfunint(tlength - t[j], gparam);
  if (t[j] <= tstart2)
  {
    gi -= gfunint(tstart2 - t[j], gparam);
  }

  double si = 0;
  for (int k = 0; k < (np - 1); ++k)
  {
    double dpx = (px[k + 1] - px[k]) / ndiv;
    double dpy = (py[k + 1] - py[k]) / ndiv;
    for (int l = 0; l < ndiv; ++l)
    {
      double x1 = px[k] + dpx * l;
      double y1 = py[k] + dpy * l;
      double x2 = px[k] + dpx * (l + 1);
      double y2 = py[k] + dpy * (l + 1);

      double det = (x1 * y2 + y1 * x[j] + x2 * y[j]) -
          (x2 * y1 + y2 * x[j] + x1 * y[j]);
      if (fabs(det) < 1.0e-10)
        continue;

      double r1 = dist(x1, y1, x[j], y[j]);
      double r2 = dist(x2, y2, x[j], y[j]);
      double phi = (r1 * r1 + r2 * r2 - dist2(x1, y1, x2, y2))/(2 * r1 * r2);
      if (fabs(phi) > 1)
        phi = 1 - 1.0e-10;

      phi = acos(phi);

      if (r1 + r2 > 1.0e-20)
      {
        double r0 = dist(x1 + r1/(r1 + r2) * (x2 - x1),
                    y1 + r1/(r1 + r2) * (y2 - y1), x[j], y[j]);

        si += sgn(det) * (ffunrint1(r1, m[j], fparam) / 6 +
            ffunrint1(r0, m[j], fparam) * 2 / 3 +
            ffunrint1(r2, m[j], fparam) / 6) * phi;
      }
    }
  }

  double intpart = kappafun(m[j], kparam) * gi * si +
      mu * integ0 / N;
  return -sumpart + intpart;
}

void etas::mloglikj1Gr(int j,
                       double mu,
                       double kparam[],
                       double gparam[],
                       double fparam[],
                       double *fvj,
                       double *dfvj)
{
  double sumpart = 0;
  double sumpartGr[8] = {0};
  if (flag[j] == 1)
  {
    double sumj = mu * bk[j];
    double sumjGr[8] = {0};
    sumjGr[0] = bk[j];

    for (int i = 0; i < j; i++)
    {
      std::array<double, 3> part1 = dkappafun(m[i], kparam);
      std::array<double, 3> part2 = dgfun(t[j] - t[i], gparam);
      std::array<double, 4> part3 = dffun1(dist2(x[j], y[j], x[i], y[i]), m[i], fparam);

      sumj    += part1[0] * part2[0] * part3[0];

      // part1_A
      sumjGr[1]  += part1[1] * part2[0] * part3[0];
      // part2_c
      sumjGr[2] += part1[0] * part2[1] * part3[0];
      //part1_alpha
      sumjGr[3]  += part1[2] * part2[0] * part3[0];
      // part2_p
      sumjGr[4] += part1[0] * part2[2] * part3[0];
      // part3_d
      sumjGr[5] += part1[0] * part2[0] * part3[1];
      // part3_q
      sumjGr[6] += part1[0] * part2[0] * part3[2];
      // part3_gamma
      sumjGr[7]  += part1[0] * part2[0] * part3[3];
    }

    sumpart = (sumj > 1.0e-25) ? log(sumj) : -100.0;

    for (int ip = 0; ip < 8; ip++)
    {
      sumpartGr[ip] += sumjGr[ip] / sumj;
    }
  }

  std::array<double, 3> int_part1 = dkappafun(m[j], kparam);

  std::array<double, 3> int_part2 = dgfunint(tlength - t[j], gparam);
  if (t[j] <= tstart2)
  {
    std::array<double, 3> gtmp = dgfunint(tstart2 - t[j], gparam);
    for (int i = 0; i < 3; i++)
    {
      int_part2[i] -= gtmp[i];
    }
  }

  double int_part3[4] = {0};
  for (int k = 0; k < (np - 1); ++k)
  {
    double dpx = (px[k + 1] - px[k]) / ndiv;
    double dpy = (py[k + 1] - py[k]) / ndiv;
    for (int l = 0; l < ndiv; ++l)
    {
      double x1 = px[k] + dpx * l;
      double y1 = py[k] + dpy * l;
      double x2 = px[k] + dpx * (l + 1);
      double y2 = py[k] + dpy * (l + 1);

      double det = (x1 * y2 + y1 * x[j] + x2 * y[j]) -
          (x2 * y1 + y2 * x[j] + x1 * y[j]);

      if (fabs(det) < 1.0e-10)
        continue;

      int id = (det < 0) ? -1 : 1;

      double r1 = dist(x1, y1, x[j], y[j]);
      double r2 = dist(x2, y2, x[j], y[j]);
      double phi = (r1 * r1 + r2 * r2 - dist2(x1, y1, x2, y2))/(2 * r1 * r2);
      if (fabs(phi) > 1)
        phi = 1 - 1.0e-10;

      phi = acos(phi);

      if (r1 + r2 > 1.0e-20)
      {
        double r0 = dist(x1 + r1/(r1 + r2) * (x2 - x1),
                    y1 + r1/(r1 + r2) * (y2 - y1), x[j], y[j]);

        std::array<double, 4> a1 = dffunrint1(r1, m[j], fparam);
        std::array<double, 4> a2 = dffunrint1(r0, m[j], fparam);
        std::array<double, 4> a3 = dffunrint1(r2, m[j], fparam);
        for (int i = 0; i < 4; i++)
          int_part3[i] += id * (a1[i] / 6 + a2[i]* 2.0 / 3 + a3[i] / 6) * phi;
      }
    }
  }

  double intpart  = int_part1[0] * int_part2[0] * int_part3[0] +
      mu * integ0 / N;
  double intpartGr[8] = {0};

  //d mu
  intpartGr[ 0 ] = integ0  / N;

  // d A
  intpartGr[ 1 ] = int_part1[1] * int_part2[0]  * int_part3[0];
  // d c
  intpartGr[ 2 ] = int_part1[0] * int_part2[1] * int_part3[0];
  // d alpha
  intpartGr[ 3 ] = int_part1[2] * int_part2[0]  * int_part3[0];
  // d p
  intpartGr[ 4 ] = int_part1[0] * int_part2[2] * int_part3[0];
  // d D
  intpartGr[ 5 ] = int_part1[0] * int_part2[0]  * int_part3[1];
  // d q
  intpartGr[ 6 ] = int_part1[0] * int_part2[0]  * int_part3[2];
  // d gamma
  intpartGr[ 7 ] = int_part1[0] * int_part2[0]  * int_part3[3];

  *fvj = -sumpart + intpart;

  for (int i = 0; i < 8; ++i)
    dfvj[i] = -sumpartGr[i] + intpartGr[i];
  return;
}


double etas::mloglikj2(int j,
                       double mu,
                       double kparam[],
                       double gparam[],
                       double fparam[])
{
  double sumpart = 0;
  if (flag[j] == 1)
  {
    double sumj = mu * bk[j];
    for (int i = 0; i < j; i++)
    {
      sumj += kappafun(m[i], kparam) *
        gfun(t[j] - t[i], gparam) *
        ffun2(dist2(x[j], y[j], x[i], y[i]), m[i], fparam);
    }

    sumpart = (sumj > 1.0e-25) ? log(sumj) : -100.0;
  }

  double gi = gfunint(tlength - t[j], gparam);
  if (t[j] <= tstart2)
  {
    gi -= gfunint(tstart2 - t[j], gparam);
  }

  double si = 0;
  for (int k = 0; k < (np - 1); ++k)
  {
    double dpx = (px[k + 1] - px[k]) / ndiv;
    double dpy = (py[k + 1] - py[k]) / ndiv;
    for (int l = 0; l < ndiv; ++l)
    {
      double x1 = px[k] + dpx * l;
      double y1 = py[k] + dpy * l;
      double x2 = px[k] + dpx * (l + 1);
      double y2 = py[k] + dpy * (l + 1);

      double det = (x1 * y2 + y1 * x[j] + x2 * y[j]) -
          (x2 * y1 + y2 * x[j] + x1 * y[j]);
      if (fabs(det) < 1.0e-10)
        continue;

      double r1 = dist(x1, y1, x[j], y[j]);
      double r2 = dist(x2, y2, x[j], y[j]);
      double phi = (r1 * r1 + r2 * r2 - dist2(x1, y1, x2, y2))/(2 * r1 * r2);
      if (fabs(phi) > 1)
        phi = 1 - 1.0e-10;

      phi = acos(phi);

      if (r1 + r2 > 1.0e-20)
      {
        double r0 = dist(x1 + r1/(r1 + r2) * (x2 - x1),
                    y1 + r1/(r1 + r2) * (y2 - y1), x[j], y[j]);

        si += sgn(det) * (ffunrint2(r1, m[j], fparam) / 6 +
            ffunrint2(r0, m[j], fparam) * 2 / 3 +
            ffunrint2(r2, m[j], fparam) / 6) * phi;
      }
    }
  }

  double intpart = kappafun(m[j], kparam) * gi * si +
      mu * integ0 / N;
  return -sumpart + intpart;
}

void etas::mloglikj2Gr(int j,
                       double mu,
                       double kparam[],
                       double gparam[],
                       double fparam[],
                       double *fvj,
                       double *dfvj)
{
  double sumpart = 0;
  double sumpartGr[7] = {0};
  if (flag[j] == 1)
  {
    double sumj = mu * bk[j];
    double sumjGr[7] = {0};
    sumjGr[0] = bk[j];

    for (int i = 0; i < j; i++)
    {
      std::array<double, 3> part1 = dkappafun(m[i], kparam);
      std::array<double, 3> part2 = dgfun(t[j] - t[i], gparam);
      std::array<double, 3> part3 = dffun2(dist2(x[j], y[j], x[i], y[i]), m[i], fparam);

      sumj    += part1[0] * part2[0] * part3[0];

      // part1_A
      sumjGr[1]  += part1[1] * part2[0] * part3[0];
      // part2_c
      sumjGr[2] += part1[0] * part2[1] * part3[0];
      //part1_alpha
      sumjGr[3]  += part1[2] * part2[0] * part3[0];
      // part2_p
      sumjGr[4] += part1[0] * part2[2] * part3[0];
      // part3_d
      sumjGr[5] += part1[0] * part2[0] * part3[1];
      // part3_gamma
      sumjGr[6]  += part1[0] * part2[0] * part3[2];
    }

    sumpart = (sumj > 1.0e-25) ? log(sumj) : -100.0;

    for (int ip = 0; ip < 7; ip++)
    {
      sumpartGr[ip] += sumjGr[ip] / sumj;
    }
  }

  std::array<double, 3> int_part1 = dkappafun(m[j], kparam);

  std::array<double, 3> int_part2 = dgfunint(tlength - t[j], gparam);
  if (t[j] <= tstart2)
  {
    std::array<double, 3> gtmp = dgfunint(tstart2 - t[j], gparam);
    for (int i = 0; i < 3; i++)
    {
      int_part2[i] -= gtmp[i];
    }
  }

  double int_part3[3] = {0};
  for (int k = 0; k < (np - 1); ++k)
  {
    double dpx = (px[k + 1] - px[k]) / ndiv;
    double dpy = (py[k + 1] - py[k]) / ndiv;
    for (int l = 0; l < ndiv; ++l)
    {
      double x1 = px[k] + dpx * l;
      double y1 = py[k] + dpy * l;
      double x2 = px[k] + dpx * (l + 1);
      double y2 = py[k] + dpy * (l + 1);

      double det = (x1 * y2 + y1 * x[j] + x2 * y[j]) -
          (x2 * y1 + y2 * x[j] + x1 * y[j]);

      if (fabs(det) < 1.0e-10)
        continue;

      int id = (det < 0) ? -1 : 1;

      double r1 = dist(x1, y1, x[j], y[j]);
      double r2 = dist(x2, y2, x[j], y[j]);
      double phi = (r1 * r1 + r2 * r2 - dist2(x1, y1, x2, y2))/(2 * r1 * r2);
      if (fabs(phi) > 1)
        phi = 1 - 1.0e-10;

      phi = acos(phi);

      if (r1 + r2 > 1.0e-20)
      {
        double r0 = dist(x1 + r1/(r1 + r2) * (x2 - x1),
                    y1 + r1/(r1 + r2) * (y2 - y1), x[j], y[j]);

        std::array<double, 3> a1 = dffunrint2(r1, m[j], fparam);
        std::array<double, 3> a2 = dffunrint2(r0, m[j], fparam);
        std::array<double, 3> a3 = dffunrint2(r2, m[j], fparam);
        for (int i = 0; i < 3; i++)
          int_part3[i] += id * (a1[i] / 6 + a2[i]* 2.0 / 3 + a3[i] / 6) * phi;
      }
    }
  }

  double intpart  = int_part1[0] * int_part2[0] * int_part3[0] +
      mu * integ0 / N;
  double intpartGr[7] = {0};

  //d mu
  intpartGr[ 0 ] = integ0  / N;

  // d A
  intpartGr[ 1 ] = int_part1[1] * int_part2[0]  * int_part3[0];
  // d c
  intpartGr[ 2 ] = int_part1[0] * int_part2[1] * int_part3[0];
  // d alpha
  intpartGr[ 3 ] = int_part1[2] * int_part2[0]  * int_part3[0];
  // d p
  intpartGr[ 4 ] = int_part1[0] * int_part2[2] * int_part3[0];
  // d D
  intpartGr[ 5 ] = int_part1[0] * int_part2[0]  * int_part3[1];
  // d gamma
  intpartGr[ 6 ] = int_part1[0] * int_part2[0]  * int_part3[2];

  *fvj = -sumpart + intpart;

  for (int i = 0; i < 7; ++i)
    dfvj[i] = -sumpartGr[i] + intpartGr[i];
  return;
}

double etas::mloglikj(int j,
                      double mu,
                      double kparam[],
                      double gparam[],
                      double fparam[])
{
  double mllj = 0;
  switch (mver)
  {
    case 1:
      mllj = mloglikj1(j, mu, kparam, gparam, fparam);
      break;
    case 2:
      mllj = mloglikj2(j, mu, kparam, gparam, fparam);
      break;
  }
  return mllj;
}

double etas::mloglik(NumericVector theta)
{
  double mu, kparam[2], gparam[2], fparam[3];
  paramhandler(theta, &mu, kparam, gparam, fparam);

  double fv = 0;

  for (int j = 0; j < N; ++j)
    fv += mloglikj(j, mu, kparam, gparam, fparam);

  return fv;
}


// ******************************************************************
// gradient of minus log likelihood function
// ******************************************************************

void etas::mloglikjGr(int j,
                      double mu,
                      double kparam[],
                      double gparam[],
                      double fparam[],
                      double *fvj,
                      double *dfvj)
{
  switch (mver)
  {
    case 1:
      mloglikj1Gr(j, mu, kparam, gparam, fparam, fvj, dfvj);
      break;
    case 2:
      mloglikj2Gr(j, mu, kparam, gparam, fparam, fvj, dfvj);
      break;
  }
}

NumericVector etas::mloglikGr(NumericVector theta)
{
  const int dimparam = theta.length();

  NumericVector out(dimparam + 1);
  double mu, kparam[2], gparam[2], fparam[3];
  paramhandler(theta, &mu, kparam, gparam, fparam);

  double fvtemp = 0, dfvtemp[8] = {};

  for (int j = 0; j < N; ++j)
  {
    double fvj, dfvj[8];
    mloglikjGr(j, mu, kparam, gparam, fparam, &fvj, dfvj);

    fvtemp += fvj;

    for (int i = 0; i < dimparam; ++i)
      dfvtemp[i] += dfvj[i];
  }

  out[0] = fvtemp;
  for (int i = 0; i < dimparam; ++i)
    out[i + 1] = dfvtemp[i] * 2 * theta[i];

  return out;
}

// ******************************************************************
// line search for the optimization algorithm
// ******************************************************************

void etas::linesearch(NumericVector xOld,
                      NumericVector h,
                      double *fv,
                      double *ram)
{
  R_CheckUserInterrupt();
  double const2 = 1.0e-16, ram1, ram2, ram3, fv1, fv2, fv3,
    a1, a2, a3, b1, b2;

  const int dimparam = xOld.length();
  
  NumericVector xNew(dimparam);
  
  if (*ram <= 1.0e-30)
    *ram = 0.1;
  
  double hnorm = 0;
  for (int i = 0; i < dimparam; i++)
    hnorm += h[i] * h[i];
  hnorm = sqrt(hnorm);

  if (hnorm > 1)
    *ram = *ram / hnorm;
  
  ram1 = 0;
  ram2 = *ram;
  fv1  = *fv;
  
  for (int i = 0; i < dimparam; i++)
    xNew[i] = xOld[i] + ram2 * h[i];
  fv2 = mloglik(xNew);

  if (fv2 > fv1)
    goto stat50;
  
  stat30:
    ram3 = ram2*2.0;
  for (int i = 0; i < dimparam ; i++)
    xNew[i] = xOld[i] + ram3 * h[i];
  fv3 = mloglik(xNew);
  if (fv3 > fv2)
    goto stat70;
  ram1 = ram2;
  ram2 = ram3;
  fv1 = fv2;
  fv2 = fv3;
  goto stat30;
  
  stat50:
    ram3 = ram2;
  fv3 = fv2;
  ram2 = ram3 * 0.1;
  if (ram2 * hnorm < const2)
  {
    *ram = 0;
    return;
  }
  for (int i = 0; i < dimparam; i++)
    xNew[i] = xOld[i] + ram2 * h[i];
  fv2 = mloglik(xNew);
  if (fv2 > fv1)
    goto stat50;
  
  stat70:
    a1 = (ram3 - ram2) * fv1;
  a2 = (ram1 - ram3) * fv2;
  a3 = (ram2 - ram1) * fv3;
  b2 = (a1 + a2 + a3) * 2;
  b1 = a1 * (ram3 + ram2) + a2 * (ram1 + ram3) + a3 * (ram2 + ram1);
  if (b2 == 0)
  {
    *ram = ram2;
    return;
  }
  else
  {
    *ram = b1 / b2;
    for (int i = 0; i < dimparam; i++)
      xNew[i] = xOld[i] + *ram*h[i];
    *fv = mloglik(xNew);
    if (*ram > ram2)
    {
      if (*fv <= fv2)
      {
        ram1 = ram2;
        ram2 = *ram;
        fv1 = fv2;
        fv2 = *fv;
        goto stat200130;
      }
      else
      {
        ram3 = *ram;
        fv3 = *fv;
        goto stat200130;
      }
    }
    else
    {
      if (*fv >= fv2)
      {
        ram1 = *ram;
        fv1 = *fv;
        goto stat200130;
      }
      else
      {
        ram3 = ram2;
        ram2 = *ram;
        fv3 = fv2;
        fv2 = *fv;
        goto stat200130;
      }
    }
  }
  
  stat200130:
    a1 = (ram3 - ram2)*fv1;
  a2 = (ram1 - ram3)*fv2;
  a3 = (ram2 - ram1)*fv3;
  b2 = (a1 + a2 + a3)*2.0;
  b1 = a1 * (ram3 + ram2) + a2 * (ram1 + ram3)
    + a3 * (ram2 + ram1);
  if (b2 == 0)
  {
    *ram = ram2;
    return;
  }
  else
  {
    *ram = b1 /b2;
    for (int i = 0; i < dimparam; i++)
      xNew[i] = xOld[i] + *ram * h[i];
    *fv = mloglik(xNew);
    if (fv2 < *fv)
      *ram = ram2;
    return;
  }
}

// ******************************************************************
// MLE for the ETAS model parameters
// ******************************************************************

List etas::fitfun(NumericVector tht,
                  NumericMatrix ihess,
                  double eps,
                  bool verbose)
{
  const int dimparam = tht.length();

  NumericVector estimate(dimparam), dfvout(dimparam);
  double fvout, aic;
  
  if (verbose)
    Rprintf("\tstart Davidon-Fletcher-Powell procedure ... \n");
  
  double tau1 = eps, tau2 = eps, eps1 = eps, eps2 = eps,
    const1 = 1.0e-17;
  
  double ramda = 0.05, fv, s1, s2;
  NumericVector s(dimparam), dx(dimparam), g0(dimparam), dg(dimparam), wrk(dimparam);
  
  // Initial estimate of inverse of hessian matrix
  NumericMatrix h = ihess;
  
  NumericVector mlkgfv = mloglikGr(tht), g(dimparam);
  fv = mlkgfv[0];
  for (int i = 0; i < dimparam; i++)
    g[i] = mlkgfv[i + 1];
  
  if (verbose)
  {
    Rprintf("Function Value = %8.4f\n", fv);
    for (int i = 0; i < dimparam; ++i)
      Rprintf("Gradient[%d] = %8.2f\ttheta[%d] = %2.6f\n", i + 1,
              g[i], i + 1, tht[i]);
  }
  
  for (int iter = 1; iter < 10; iter++)
  {
    R_CheckUserInterrupt();
    for (int ic = 0; ic < dimparam; ic++)
    {
      if (ic > 0 || iter > 1)
      {
        for (int i = 0; i < dimparam; i++)
          dg[i] = g[i] - g0[i];
        
        for (int i = 0; i < dimparam; i++)
        {
          double sum = 0;
          for (int j = 0; j < dimparam; j++)
            sum += dg[j] * h(i, j);
          wrk[i] = sum;
        }
        
        s1 = 0.0;
        s2 = 0.0;
        for (int i = 0; i < dimparam; i++)
        {
          s1 += wrk[i] * dg[i];
          s2 += dx[i] * dg[i];
        }
        
        if (s1 <= const1 || s2 <= const1)
        {
          fvout = -fv;
          aic = 2 * (fv + dimparam);
          if (verbose)
            Rprintf ("loglikelihood = %8.5f\tAIC = %8.5f\n", -fv, aic);
          for( int i = 0; i < dimparam; i++ )
          {
            dfvout[i] = g[i];
            estimate[i] = tht[i];
            for (int j = 0; j < dimparam; j++)
              ihess(i, j) = h(i, j);
            if (verbose)
              Rprintf("theta[%d] = %2.8f\t gradient[%d] = %8.4f\n",
                      i + 1, pow(tht[i], 2), i + 1, g[i]);
          }
          
          return List::create(Named("estimate") = estimate,
                              Named("fvout") = fvout,
                              Named("dfvout") = dfvout,
                              Named("aic") = aic,
                              Named("hess") = ihess);
        }
        
        if (s1 <= s2)
        {
          // fletcher type correction
          for (int i = 0; i < dimparam; i++)
            for (int j = i; j < dimparam; j++)
            {
              h(i, j) -= (dx[i] * wrk[j] + wrk[i] * dx[j] -
                dx[i] * dx[j] * (1 + s1 / s2)) / s2;
              h(j, i) = h(i, j);
            }
        }
        else
        {
          // Update the inverse of hessian matrix
          for (int i = 0; i < dimparam; i++)
            for (int j = i; j < dimparam; j++)
            {	// davidon-fletcher-powell type correction
              h(i, j) += dx[i] * dx[j] / s2 - wrk[i] * wrk[j] / s1;
              h(j, i) = h(i, j);
            }
        }
      }
      
      double ss = 0;
      for (int i = 0; i < dimparam; i++)
      {
        double sum = 0;
        for (int j = 0; j < dimparam; j++)
          sum += h(i, j) * g[j];
        ss += sum * sum;
        s[i] = -sum;
      }
      s1 = 0.0;
      s2 = 0.0;
      for (int i = 0; i < dimparam; i++)
      {
        s1 += s[i] * g[i];
        s2 += g[i] * g[i];
      }
      
      if ((fabs(s1) / sqrt(s2) <= tau1) && (sqrt(s2) <= tau2))
      {
        fvout = -fv;
        aic = 2 * (fv + dimparam);
        if (verbose)
          Rprintf ("loglikelihood = %8.5f\tAIC = %8.5f\n", -fv, aic);
        for( int i = 0; i < dimparam; i++ )
        {
          dfvout[i] = g[i];
          estimate[i] = tht[i];
          for (int j = 0; j < dimparam; j++)
            ihess(i, j) = h(i, j);
          if (verbose)
            Rprintf("theta[%d] = %2.8f\t gradient[%d] = %8.4f\n",
                    i + 1, pow(tht[i], 2), i + 1, g[i]);
        }
        
        return List::create(Named("estimate") = estimate,
                            Named("fvout") = fvout,
                            Named("dfvout") = dfvout,
                            Named("aic") = aic,
                            Named("hess") = ihess);
      }
      
      if (s1 >= 0)
      {
        for (int i = 0; i < dimparam; i++)
        {
          for (int j = 0; j < dimparam; j++)
            h(i, j) = 0.0;
          h(i, i) = 1.0;
          s[i] = -s[i];
        }
      }
      
      double ed = fv;
      if (verbose)
        Rprintf("\nline search along the specified direction ...");
      // line  search
      linesearch(tht, s, &ed, &ramda);
      
      if (verbose)
        Rprintf(" zeta = %f\n", ramda);
      
      //R_CheckUserInterrupt();
      
      s1 = 0;
      for (int i = 0; i < dimparam; i++)
      {
        dx[i] = s[i] * ramda;
        s1 += dx[i] * dx[i];
        g0[i] = g[i];
        tht[i] += dx[i];
      }
      
      double fv0 = fv;
      mlkgfv = mloglikGr(tht);
      fv = mlkgfv[0];
      for (int i = 0; i < dimparam; i++)
        g[i] = mlkgfv[i + 1];
      
      if (verbose)
      {
        Rprintf("Function Value = %8.4f\n", fv);
        for (int i = 0; i < dimparam; ++i)
          Rprintf("Gradient[%d] = %8.2f\ttheta[%d] = %2.6f\n", i + 1,
                  g[i], i + 1, tht[i]);
      }
      
      s2 = 0;
      for (int i = 0; i < dimparam; i++)
        s2 += g[i] * g[i];
      if (sqrt(s2) > tau2)
        continue;
      
      if (fv0/fv - 1 < eps1 && sqrt(s1) < eps2)
      {
        fvout = -fv;
        aic = 2 * (fv + dimparam);
        if (verbose)
          Rprintf ("loglikelihood = %8.5f\tAIC = %8.5f\n", -fv, aic);
        for( int i = 0; i < dimparam; i++ )
        {
          dfvout[i] = g[i];
          estimate[i] = tht[i];
          for (int j = 0; j < dimparam; j++)
            ihess(i, j) = h(i, j);
          if (verbose)
            Rprintf("theta[%d] = %2.8f\t gradient[%d] = %8.4f\n", i + 1, pow(tht[i], 2), i + 1, g[i]);
        }
        
        return List::create(Named("estimate") = estimate,
                            Named("fvout") = fvout,
                            Named("dfvout") = dfvout,
                            Named("aic") = aic,
                            Named("hess") = ihess);
      }
    }
  }
  return 0;
}

// ******************************************************************
// minus log likelihood function: parallel computing
// ******************************************************************

double etas::mloglikMP(NumericVector theta,
                       int nthreads)
{
  double mu, kparam[2], gparam[2], fparam[3];
  paramhandler(theta, &mu, kparam, gparam, fparam);

  double fv = 0;
  
  #pragma omp parallel num_threads(nthreads)
  {
    double fv_thread = 0;
  
    #pragma omp for
    for (int j = 0; j < N; ++j)
    {
      fv_thread += mloglikj(j, mu, kparam, gparam, fparam);
    }

    #pragma omp critical
    {
      fv += fv_thread;
    }
  }

  return fv;
}

// ******************************************************************
// gradient of minus log likelihood function: parallel computing
// ******************************************************************

NumericVector etas::mloglikGrMP(NumericVector theta, int nthreads)
{
  const int dimparam = theta.length();

  NumericVector out(dimparam + 1);

  double mu, kparam[2], gparam[2], fparam[3];
  paramhandler(theta, &mu, kparam, gparam, fparam);

  double fvtemp = 0, dfvtemp[8] = {};
  
  #pragma omp parallel num_threads(nthreads)
  {
    double fvtemp_thread = 0, dfvtemp_thread[8] = {};
  
    #pragma omp for //schedule(static)
    for (int j = 0; j < N; ++j)
    {
      double fvj, dfvj[8];
      mloglikjGr(j, mu, kparam, gparam, fparam, &fvj, dfvj);

      fvtemp_thread += fvj;
      for (int i = 0; i < dimparam; ++i)
        dfvtemp_thread[i] += dfvj[i];
    }

    #pragma omp critical
    {
      fvtemp += fvtemp_thread;
      for (int i = 0; i < dimparam; ++i)
        dfvtemp[i] += dfvtemp_thread[i];
    }
  }

  out[0] = fvtemp;
  for (int i = 0; i < dimparam; ++i)
    out[i + 1] = dfvtemp[i]  * 2 * theta[i];

  return out;
}


// ******************************************************************
// line search for the optimization algorithm: parallel computing
// ******************************************************************

void etas::linesearchMP(NumericVector xOld,
                        NumericVector h,
                        double *fv,
                        double *ram,
                        int nthreads)
{
  R_CheckUserInterrupt();
  double const2 = 1.0e-16, ram1, ram2, ram3, fv1, fv2, fv3,
    a1, a2, a3, b1, b2;

  const int dimparam = xOld.length();

  NumericVector xNew(dimparam);
  
  if (*ram <= 1.0e-30)
    *ram = 0.1;
  
  double hnorm = 0;
  for (int i = 0; i < dimparam; i++)
    hnorm += h[i] * h[i];
  hnorm = sqrt(hnorm);

  if (hnorm > 1)
    *ram = *ram / hnorm;
  
  ram1 = 0;
  ram2 = *ram;
  fv1  = *fv;
  
  for (int i = 0; i < dimparam; i++)
    xNew[i] = xOld[i] + ram2 * h[i];
  fv2 = mloglikMP(xNew, nthreads);
  
  if (fv2 > fv1)
    goto stat50;
  
  stat30:
    ram3 = ram2*2.0;
  for (int i = 0; i < dimparam; i++)
    xNew[i] = xOld[i] + ram3 * h[i];
  fv3 = mloglikMP(xNew, nthreads);
  if (fv3 > fv2)
    goto stat70;
  ram1 = ram2;
  ram2 = ram3;
  fv1 = fv2;
  fv2 = fv3;
  goto stat30;
  
  stat50:
    ram3 = ram2;
  fv3 = fv2;
  ram2 = ram3 * 0.1;
  if (ram2 * hnorm < const2)
  {
    *ram = 0;
    return;
  }
  for (int i = 0; i < dimparam; i++)
    xNew[i] = xOld[i] + ram2 * h[i];
  fv2 = mloglikMP(xNew, nthreads);
  if (fv2 > fv1)
    goto stat50;
  
  stat70:
    a1 = (ram3 - ram2) * fv1;
  a2 = (ram1 - ram3) * fv2;
  a3 = (ram2 - ram1) * fv3;
  b2 = (a1 + a2 + a3) * 2;
  b1 = a1 * (ram3 + ram2) + a2 * (ram1 + ram3) + a3 * (ram2 + ram1);
  if (b2 == 0)
  {
    *ram = ram2;
    return;
  }
  else
  {
    *ram = b1 / b2;
    for (int i = 0; i < dimparam; i++)
      xNew[i] = xOld[i] + *ram*h[i];
    *fv = mloglikMP(xNew, nthreads);
    if (*ram > ram2)
    {
      if (*fv <= fv2)
      {
        ram1 = ram2;
        ram2 = *ram;
        fv1 = fv2;
        fv2 = *fv;
        goto stat200130;
      }
      else
      {
        ram3 = *ram;
        fv3 = *fv;
        goto stat200130;
      }
    }
    else
    {
      if (*fv >= fv2)
      {
        ram1 = *ram;
        fv1 = *fv;
        goto stat200130;
      }
      else
      {
        ram3 = ram2;
        ram2 = *ram;
        fv3 = fv2;
        fv2 = *fv;
        goto stat200130;
      }
    }
  }
  
  stat200130:
    a1 = (ram3 - ram2)*fv1;
  a2 = (ram1 - ram3)*fv2;
  a3 = (ram2 - ram1)*fv3;
  b2 = (a1 + a2 + a3)*2.0;
  b1 = a1 * (ram3 + ram2) + a2 * (ram1 + ram3)
    + a3 * (ram2 + ram1);
  if (b2 == 0)
  {
    *ram = ram2;
    return;
  }
  else
  {
    *ram = b1 /b2;
    for (int i = 0; i < dimparam; i++)
      xNew[i] = xOld[i] + *ram*h[i];
    *fv = mloglikMP(xNew, nthreads);
    if (fv2 < *fv)
      *ram = ram2;
    return;
  }
}

// ******************************************************************
// MLE for the ETAS model parameters
// ******************************************************************

List etas::fitfunMP(NumericVector tht,
                    NumericMatrix ihess,
                    double eps,
                    bool verbose,
                    int nthreads)
{
  const int dimparam = tht.length();

  NumericVector estimate(dimparam), dfvout(dimparam);
  double fvout, aic;
  
  if (verbose)
    Rprintf("\tstart Davidon-Fletcher-Powell procedure ... \n");
  
  double tau1 = eps, tau2 = eps, eps1 = eps, eps2 = eps, const1 = 1.0e-17;
  
  double ramda = 0.05, fv, s1, s2;
  NumericVector s(dimparam), dx(dimparam), g0(dimparam), dg(dimparam), wrk(dimparam);

  // Initial estimate of inverse of hessian matrix
  NumericMatrix h = ihess;
  
  NumericVector mlkgfv = mloglikGrMP(tht, nthreads), g(dimparam);
  fv = mlkgfv[0];
  for (int i = 0; i < dimparam; i++)
    g[i] = mlkgfv[i + 1];
  
  if (verbose)
  {
    Rprintf("Function Value = %8.4f\n", fv);
    for (int i = 0; i < dimparam; ++i)
      Rprintf("Gradient[%d] = %8.2f\ttheta[%d] = %2.6f\n",
              i + 1, g[i], i + 1, tht[i]);
  }
  
  for (int iter = 1; iter < 10; iter++)
  {
    R_CheckUserInterrupt();
    for (int ic = 0; ic < dimparam; ic++)
    {
      if (ic > 0 || iter > 1)
      {
        for (int i = 0; i < dimparam; i++)
          dg[i] = g[i] - g0[i];
        
        for (int i = 0; i < dimparam; i++)
        {
          double sum = 0;
          for (int j = 0; j < dimparam; j++)
            sum += dg[j] * h(i, j);
          wrk[i] = sum;
        }
        
        s1 = 0.0;
        s2 = 0.0;
        for (int i = 0; i < dimparam; i++)
        {
          s1 += wrk[i] * dg[i];
          s2 += dx[i] * dg[i];
        }
        
        if (s1 <= const1 || s2 <= const1)
        {
          fvout = -fv;
          aic = 2 * (fv + dimparam);
          if (verbose)
            Rprintf ("loglikelihood = %8.5f\tAIC = %8.5f\n", -fv, aic);
          for( int i = 0; i < dimparam; i++ )
          {
            dfvout[i] = g[i];
            estimate[i] = tht[i];
            for (int j = 0; j < dimparam; j++)
              ihess(i, j) = h(i, j);
            if (verbose)
              Rprintf("theta[%d] = %2.8f\t gradient[%d] = %8.4f\n",
                      i + 1, pow(tht[i], 2), i + 1, g[i]);
          }
          
          return List::create(Named("estimate") = estimate,
                              Named("fvout") = fvout,
                              Named("dfvout") = dfvout,
                              Named("aic") = aic,
                              Named("hess") = ihess);
        }
        
        if (s1 <= s2)
        {
          // fletcher type correction
          for (int i = 0; i < dimparam; i++)
            for (int j = i; j < dimparam; j++)
            {
              h(i, j) -= (dx[i] * wrk[j] + wrk[i] * dx[j] -
                dx[i] * dx[j] * (1 + s1 / s2)) / s2;
              h(j, i) = h(i, j);
            }
        }
        else
        {
          // Update the inverse of hessian matrix
          for (int i = 0; i < dimparam; i++)
            for (int j = i; j < dimparam; j++)
            {	// davidon-fletcher-powell type correction
              h(i, j) += dx[i] * dx[j]/s2 - wrk[i] * wrk[j] / s1;
              h(j, i) = h(i, j);
            }
        }
      }
      double ss = 0;
      for (int i = 0; i < dimparam; i++)
      {
        double sum = 0;
        for (int j = 0; j < dimparam; j++)
          sum += h(i, j) * g[j];
        ss += sum * sum;
        s[i] = -sum;
      }
      s1 = 0.0;
      s2 = 0.0;
      for (int i = 0; i < dimparam; i++)
      {
        s1 += s[i] * g[i];
        s2 += g[i] * g[i];
      }
      
      if ((fabs(s1) / sqrt(s2) <= tau1) && (sqrt(s2) <= tau2))
      {
        fvout = -fv;
        aic = 2 * (fv + dimparam);
        if (verbose)
          Rprintf ("loglikelihood = %8.5f\tAIC = %8.5f\n",
                   -fv, aic);
        for( int i = 0; i < dimparam; i++ )
        {
          dfvout[i] = g[i];
          estimate[i] = tht[i];
          for (int j = 0; j < dimparam; j++)
            ihess(i, j) = h(i, j);
          if (verbose)
            Rprintf("theta[%d] = %2.8f\t gradient[%d] = %8.4f\n",
                    i + 1, pow(tht[i], 2), i + 1, g[i]);
        }
        
        return List::create(Named("estimate") = estimate,
                            Named("fvout") = fvout,
                            Named("dfvout") = dfvout,
                            Named("aic") = aic,
                            Named("hess") = ihess);
      }
      
      if (s1 >= 0)
        for (int i = 0; i < dimparam; i++)
        {
          for (int j = 0; j < dimparam; j++)
            h(i, j) = 0.0;
          h(i, i) = 1.0;
          s[i] = -s[i];
        }
        
      double ed = fv;
      if (verbose)
        Rprintf("\nline search along the specified direction ...");
      // line  search
      linesearchMP(tht, s, &ed, &ramda, nthreads);
      
      if (verbose)
        Rprintf(" zeta = %f\n", ramda);
      
      //R_CheckUserInterrupt();
      
      s1 = 0;
      for (int i = 0; i < dimparam; i++)
      {
        dx[i] = s[i] * ramda;
        s1 += dx[i] * dx[i];
        g0[i] = g[i];
        tht[i] += dx[i];
      }
      
      double fv0 = fv;
      mlkgfv = mloglikGrMP(tht, nthreads);
      fv = mlkgfv[0];
      for (int i = 0; i < dimparam; i++)
        g[i] = mlkgfv[i + 1];
      
      if (verbose)
      {
        Rprintf("Function Value = %8.4f\n", fv);
        for (int i = 0; i < dimparam; ++i)
          Rprintf("Gradient[%d] = %8.2f\ttheta[%d] = %2.6f\n",
                  i + 1, g[i], i + 1, tht[i]);
      }
      
      s2 = 0;
      for (int i = 0; i < dimparam; i++)
        s2 += g[i] * g[i];
      if (sqrt(s2) > tau2)
        continue;
      if (fv0/fv - 1 < eps1 && sqrt(s1) < eps2)
      {
        fvout = -fv;
        aic = 2 * (fv + dimparam);
        if (verbose)
          Rprintf ("loglikelihood = %8.5f\tAIC = %8.5f\n", -fv, aic);
        for( int i = 0; i < dimparam; i++ )
        {
          dfvout[i] = g[i];
          estimate[i] = tht[i];
          for (int j = 0; j < dimparam; j++)
            ihess(i, j) = h(i, j);
          if (verbose)
            Rprintf("theta[%d] = %2.8f\t gradient[%d] = %8.4f\n", i + 1, pow(tht[i], 2), i + 1, g[i]);
        }
        
        return List::create(Named("estimate") = estimate,
                            Named("fvout") = fvout,
                            Named("dfvout") = dfvout,
                            Named("aic") = aic,
                            Named("hess") = ihess);
      }
    }
  }
  return 0;
}


// ******************************************************************
// wrapper fit function for R
// ******************************************************************

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
List cxxfit(NumericVector tht,
            NumericMatrix revents,
            NumericMatrix rpoly,
            NumericVector tperiod,
            double rinteg0,
            NumericMatrix ihess,
            int ndiv,
            double eps,
            bool verbose,
            int nthreads,
            int mver)
{
  etas data;
  data.set(revents, rpoly, tperiod, rinteg0, ndiv, mver);
  
  #ifdef _OPENMP
  if (nthreads > 1)
  {
    //setenv("OMP_STACKSIZE", "200M", 1);
    omp_set_dynamic(0);
    return data.fitfunMP(tht, ihess, eps, verbose, nthreads);
  }
  else
    return data.fitfun(tht, ihess, eps, verbose);
  #else
  // serial version of code
  return data.fitfun(tht, ihess, eps, verbose);
  #endif
}

// *******************************************************************************
// *******************************************************************************
// *******************************************************************************
// *******************************************************************************
// *******************************************************************************
// *******************************************************************************
// *******************************************************************************
// *******************************************************************************
// *******************************************************************************

// ******************************************************************
// conditional intensity function at (tv xv, yv)
// ******************************************************************
NumericVector lambda(NumericVector tv,
                     NumericVector xv,
                     NumericVector yv,
                     NumericVector theta,
                     NumericMatrix revents)
{
  NumericVector t = revents( _, 0), x = revents( _, 1), y = revents( _, 2),
    m = revents( _, 3), bk = revents( _, 5);
  
  //const double mu = theta[0];
  const double A = theta[1];
  const double c = theta[2];
  const double alpha = theta[3];
  const double p = theta[4];
  const double D =  theta[5];
  const double q=  theta[6];
  const double gamma = theta[7];

  double kparam[] = {A, alpha};
  double gparam[] = {c, p};
  double fparam[] = {D, gamma, q};

  NumericVector out(tv.length());
  double s = 0; //= mu * bk[j];
  
  for (int j = 0; j < tv.length(); j++)
  {
    int i = 0;
    while (t[i] < tv[j])
    {
      s += kappafun(m[i], kparam) *
        gfun(tv[j] - t[i], gparam) *
        ffun1(dist2(xv[j], yv[j], x[i], y[i]), m[i], fparam);
      i++;
    }
    out[j] = s;
  }
  
  return out;
}
// ******************************************************************

// *******************************************************************************
// *******************************************************************************
// *******************************************************************************
// *******************************************************************************
// *******************************************************************************
// *******************************************************************************
// *******************************************************************************
// *******************************************************************************
// *******************************************************************************
// *******************************************************************************


// ******************************************************************
// bivariate Gaussian kernel function with bandwidth sig * I_2
// ******************************************************************

inline
  double dGauss(double r2, double sig)
  {
    return exp(-r2 /(2 * sig * sig)) / (2 * M_PI * sig * sig);
  }

// integral of the gaussian kernel with bandwidth w[0] from 0 to r
inline
  double pGauss(double r, double w)
  {
    return (1 - exp(-(r * r) / (2 * w * w))) / (2 * M_PI);
  }

// ******************************************************************
// stochastic declustring algorithm
// ******************************************************************
// [[Rcpp::export]]
List cxxdeclust(NumericVector param,
                NumericMatrix revents,
                NumericMatrix rpoly,
                NumericVector bwd,
                NumericVector tperiod,
                int ndiv,
                int mver)
{
  NumericVector t = revents( _, 0), x = revents( _, 1), y = revents( _, 2),
    m = revents( _, 3), flag = revents( _, 4), bk = revents( _, 5),
    pb = revents( _, 6), lam = revents( _, 7);
  NumericVector px = rpoly( _, 0), py = rpoly( _, 1);
  
  int N = t.length(), np = px.length();
  
  // extract time period information
  const double tstart2 = tperiod[0], tlength = tperiod[1];

  modelhandler model;
  model.set(mver, param);

  double integ0 = 0;
  for (int i = 0; i < N; i++)
  {
    double s = 0;
    for (int j = 0; j < N; j++)
    {
      s += pb[j] * dGauss(dist2(x[i], y[i], x[j], y[j]), bwd[j]);
    }
    bk[i] = s / (tlength - tstart2);
    
    double sum = 0, dpx, dpy, x1, x2, y1, y2, det, r0, r1, r2, phi;
    
    for (int k = 0; k < (np - 1); ++k)
    {
      dpx = (px[k + 1] - px[k]) / ndiv;
      dpy = (py[k + 1] - py[k]) / ndiv;
      for (int l = 0; l < ndiv; ++l)
      {
        x1 = px[k] + dpx * l;
        y1 = py[k] + dpy * l;
        x2 = px[k] + dpx * (l + 1);
        y2 = py[k] + dpy * (l + 1);
        
        det = (x1 * y2 + y1 * x[i] + x2 * y[i]) -
          (x2 * y1 + y2 * x[i] + x1 * y[i]);
        
        if (fabs(det) < 1.0e-10)
          continue;
        
        r1 = dist(x1, y1, x[i], y[i]);
        r2 = dist(x2, y2, x[i], y[i]);
        phi = (r1 * r1 + r2 * r2 - dist2(x1, y1, x2, y2))/(2 * r1 * r2);
        if (fabs(phi) > 1)
          phi = 1 - 1.0e-10;
        
        phi = acos(phi);
        
        if (r1 + r2 > 1.0e-20)
        {
          r0 = dist(x1 + r1/(r1 + r2) * (x2 - x1),
                    y1 + r1/(r1 + r2) * (y2 - y1), x[i], y[i]);
          
          sum += sgn(det) * (pGauss(r1, bwd[i])/6 + (pGauss(r0, bwd[i]) * 2)/3 +
            pGauss(r2, bwd[i])/6) * phi;
        }
      }
    }
    
    integ0 += pb[i] * sum;
    
    double s_thread = model.mufun() * bk[i];
    for (int j = 0; j < i; ++j)
    {
      s_thread += model.kappafun0(m[j]) *
        model.gfun0(t[i] - t[j]) *
        model.ffun0(dist2(x[i], y[i], x[j], y[j]), m[j]);
    }
    
    lam[i] = s_thread;
    
    revents(i, 5) = bk[i];
    revents(i, 6) = model.mufun() * bk[i] / lam[i]; // probability of event i being a background event
    revents(i, 7) = lam[i];
  }
  
  return List::create(Named("revents") = revents,
                      Named("integ0") = integ0);
}

// ******************************************************************
// output rates algorithm
// ******************************************************************
// [[Rcpp::export]]
List cxxrates(NumericVector param,
              NumericMatrix revents,
              NumericVector bwd,
              NumericVector tperiod,
              NumericVector gx,
              NumericVector gy,
              int mver)
{
  NumericVector t = revents( _, 0), x = revents( _, 1), y = revents( _, 2),
    m = revents( _, 3), pb = revents( _, 6);
  
  // extract time period information
  const double tstart2 = tperiod[0], tlength = tperiod[1];

  modelhandler model;
  model.set(mver, param);

  int N = t.length(), ngx = gx.length(), ngy = gy.length();
  
  NumericMatrix bkgd(ngx, ngy), total(ngx, ngy), clust(ngx, ngy),
  lamb(ngx, ngy);
  
  double tmp, sum1, sum2;
  for (int i = 0; i < ngx; i++)
    for (int j = 0; j < ngy; j++)
    {
      sum1 = sum2 = 0;
      for (int l = 0; l < N; l++)
      {
        tmp = dGauss(dist2(x[l], y[l], gx[i], gy[j]), bwd[l]);
        sum1 += pb[l] * tmp;
        sum2 += tmp;
      }
      bkgd(i, j) = sum1 / (tlength - tstart2);
      total(i, j) = sum2 / (tlength - tstart2);
      clust(i, j) = 1 - sum1 / sum2;
      lamb(i, j) = model.mufun() * bkgd(i, j);
      
      for (int l = 0; l < N; l++)
      {
        lamb(i, j) += model.kappafun0(m[l]) *
          model.gfun0(tlength - t[l]) *
          model.ffun0(dist2(x[l], y[l], gx[i], gy[j]), m[l]);
      }
    }
  
  return List::create(Named("bkgd") = bkgd,
                      Named("total") = total,
                      Named("clust") = clust,
                      Named("lamb") = lamb);
}

// ******************************************************************
// transformed times: \tau_i
// ******************************************************************

// [[Rcpp::export]]
NumericVector cxxtimetrans(NumericVector theta,
                           NumericMatrix revents,
                           NumericMatrix rpoly,
                           NumericVector tperiod,
                           double integ0,
                           int ndiv,
                           int mver)
{
  NumericVector t = revents( _, 0), x = revents( _, 1), y = revents( _, 2),
    m = revents( _, 3);
  NumericVector px = rpoly( _, 0), py = rpoly( _, 1);
  const double tstart2 = tperiod[0], tlength = tperiod[1];

  modelhandler model;
  model.set(mver, theta);

  const int N = revents.nrow();
  NumericVector sinteg(N), out(N);
  
  for (int i=0; i < N; i++)
  {
    double si = 0;
    for (int k = 0; k < (px.length() - 1); ++k)
    {
      double dxx = (px[k + 1] - px[k]) / ndiv;
      double dyy = (py[k + 1] - py[k]) / ndiv;
      for (int l = 0; l < ndiv; ++l)
      {
        double x1 = px[k] + dxx * l;
        double y1 = py[k] + dyy * l;
        double x2 = px[k] + dxx * (l + 1);
        double y2 = py[k] + dyy * (l + 1);
        double det = (x1 * y2 + y1 * x[i] + x2 * y[i]) -
        (x2 * y1 + y2 * x[i] + x1 * y[i]);

        if (fabs(det) < 1.0e-10)
          continue;

        double r1 = dist(x1, y1, x[i], y[i]);
        double r2 = dist(x2, y2, x[i], y[i]);
        double phi = (r1 * r1 + r2 * r2 - dist2(x1, y1, x2, y2))/(2 * r1 * r2);
        if (fabs(phi) > 1)
          phi = 1 - 1.0e-10;

        phi = acos(phi);

        if (r1 + r2 > 1.0e-20)
        {
          double r0 = dist(x1 + r1/(r1 + r2) * (x2 - x1),
                           y1 + r1/(r1 + r2) * (y2 - y1), x[i], y[i]);

          si += sgn(det) * (model.ffunrint0(r1, m[i]) / 6 +
            model.ffunrint0(r0, m[i]) * 2 / 3 +
            model.ffunrint0(r2, m[i]) / 6) * phi;
        }
      }
    }

    sinteg[i] = model.kappafun0(m[i]) * si;
  }
  
  for (int j=0; j < N; ++j)
  {
    double sum = 0;
    for (int i=0; i < j; i++)
    {
      if (t[i] > tstart2)
        sum += model.gfunint0(t[j] - t[i]) * sinteg[i];
      else
        sum += (model.gfunint0(t[j] - t[i]) -
          model.gfunint0(tstart2 - t[i])) * sinteg[i];
    }
    out[j] = model.mufun() * integ0 * (t[j] - tstart2) / (tlength - tstart2) + sum;
  }
  return out;
}

// ******************************************************************
// temporal intensity function: integrating over the spatial domain
// ******************************************************************

// [[Rcpp::export]]
NumericVector cxxlambdtemp(NumericVector tg,
                           NumericVector theta,
                           NumericMatrix revents,
                           NumericMatrix rpoly,
                           NumericVector tperiod,
                           double integ0,
                           int ndiv,
                           int mver)
{
  NumericVector t = revents( _, 0), x = revents( _, 1), y = revents( _, 2),
    m = revents( _, 3);
  NumericVector px = rpoly( _, 0), py = rpoly( _, 1);
  const double tstart2 = tperiod[0], tlength = tperiod[1];
  
  modelhandler model;
  model.set(mver, theta);

  const int N = revents.nrow();
  NumericVector sinteg(N);
  const int ng = tg.length();
  NumericVector out(ng);
  
  for (int i=0; i < N; i++)
  {
    double si = 0;
    for (int k = 0; k < (px.length() - 1); ++k)
    {
      double dxx = (px[k + 1] - px[k]) / ndiv;
      double dyy = (py[k + 1] - py[k]) / ndiv;
      for (int l = 0; l < ndiv; ++l)
      {
        double x1 = px[k] + dxx * l;
        double y1 = py[k] + dyy * l;
        double x2 = px[k] + dxx * (l + 1);
        double y2 = py[k] + dyy * (l + 1);
        double det = (x1 * y2 + y1 * x[i] + x2 * y[i]) -
        (x2 * y1 + y2 * x[i] + x1 * y[i]);

        if (fabs(det) < 1.0e-10)
          continue;

        double r1 = dist(x1, y1, x[i], y[i]);
        double r2 = dist(x2, y2, x[i], y[i]);
        double phi = (r1 * r1 + r2 * r2 - dist2(x1, y1, x2, y2))/(2 * r1 * r2);
        if (fabs(phi) > 1)
          phi = 1 - 1.0e-10;

        phi = acos(phi);

        if (r1 + r2 > 1.0e-20)
        {
          double r0 = dist(x1 + r1/(r1 + r2) * (x2 - x1),
                           y1 + r1/(r1 + r2) * (y2 - y1), x[i], y[i]);

          si += sgn(det) * (model.ffunrint0(r1, m[i]) / 6 +
            model.ffunrint0(r0, m[i]) * 2 / 3 +
            model.ffunrint0(r2, m[i]) / 6) * phi;
        }
      }
    }
    sinteg[i] = model.kappafun0(m[i]) * si;
  }
  
  for (int j=0; j < ng; ++j)
  {
    double sum = 0;
    for (int i=0; i < N; i++)
    {
      if (t[i] < tg[j])
      {
        sum += model.gfun0(tg[j] - t[i]) * sinteg[i];
      }
    }
    out[j] = model.mufun() * integ0 /(tlength - tstart2) + sum;
  }
  return out;
}

// ******************************************************************
// spatial intensity function: integrating over the temporal domain
// ******************************************************************

// [[Rcpp::export]]
NumericVector cxxlambspat(NumericVector xg,
                          NumericVector yg,
                          NumericVector theta,
                          NumericMatrix revents,
                          NumericMatrix rpoly,
                          NumericVector tperiod,
                          NumericVector bwd,
                          int mver)
{
  NumericVector t = revents( _, 0), x = revents( _, 1), y = revents( _, 2),
    m = revents( _, 3), bk = revents( _, 5), pb = revents( _, 6);
  NumericVector px = rpoly( _, 0), py = rpoly( _, 1);
  const double tstart2 = tperiod[0], tlength = tperiod[1];

  modelhandler model;
  model.set(mver, theta);

  const int N = revents.nrow();
  
  const int ng = xg.length();
  NumericVector out(ng);
  
  for (int j=0; j < ng; ++j)
  {
    double sum = 0, s1 = 0, s2 = 0, gint;
    for (int i = 0; i < N; i++)
    {
      if (t[i] > tstart2)
      {
        gint = model.gfunint0(tlength - t[i]);
      }
      else
      {
        gint = model.gfunint0(tlength - t[i]) -
          model.gfunint0(tstart2 - t[i]);
      }
      double r2 = dist2(xg[j], yg[j], x[i], y[i]);
      sum += model.kappafun0(m[i]) * gint * model.ffun0(r2, m[i]);
      s1 += exp(-r2/(2 * bwd[i] * bwd[i])) / (2 * M_PI * bwd[i] * bwd[i]);
      s2 += pb[i] *  s1;
    }
    
    out[j] =  sum + model.mufun() * s2/(tlength - tstart2);
  }
  
  return out;
}

// ******************************************************************
// variable bandwidth kernel smoothing
// ******************************************************************
// [[Rcpp::export]]
List cxxSmooth(NumericVector x,
               NumericVector y,
               NumericVector bwd,
               NumericVector gx,
               NumericVector gy,
               bool expand)
{

  const int N = x.length(), ngx = gx.length(), ngy = gy.length();
  double sum;
  
  if (expand)
  {
    NumericMatrix out(ngx, ngy);
    for (int i = 0; i < ngx; i++)
    {
      R_CheckUserInterrupt();
      for (int j = 0; j < ngy; j++)
      {
        sum = 0;
        for (int l = 0; l < N; l++)
        {
          sum += dGauss(dist2(x[l], y[l], gx[i], gy[j]), bwd[l]);
        }
        out(i, j) = sum;
      }
    }
    return List::create(Named("out") = out);
  }
  else{
    if (ngx != ngy)
      stop("gird coordinates must have the same length.");
    NumericVector out(ngx);
    for (int i = 0; i < ngx; i++)
    {
      sum = 0;
      for (int l = 0; l < N; l++)
      {
        sum += dGauss(dist2(x[l], y[l], gx[i], gy[i]), bwd[l]);
      }
      out[i] = sum;
    }
    return List::create(Named("out") = out);
  }
}

// ******************************************************************
// distant based test statistic for the spatio-temporal Poisson test
// ******************************************************************
// [[Rcpp::export]]
double cxxstpoisstest(NumericVector xrank, 
                      NumericVector yrank,
                      NumericMatrix M)
{
  const int n = xrank.length();
  NumericMatrix tmp(n, n);
  double dfv, dfvtmp, out = 0;
  
  for (int k = 1; k <= n; k++)
  {
    R_CheckUserInterrupt();
    dfv = 0;
    for (int i = 0; i < n; i++)
    {
      for (int j =0; j < n; j++)
      {
        if ((yrank[i] >= yrank[k - 1]) && (xrank[j] >= xrank[k - 1]))
          tmp(i, j) += 1;
        dfvtmp = tmp(i, j)/n - M(i, j)/n * k/n;
        if (fabs(dfvtmp) > dfv)
          dfv = dfvtmp;
      }
    }
    if (dfv > out)
      out = dfv;
  }
  return out;
}

// Parallel version
// [[Rcpp::export]]
double cxxstpoisstestMP(NumericVector xrank, 
                      NumericVector yrank,
                      NumericMatrix M,
                      int nthreads)
{
  const int n = xrank.length();
  NumericMatrix tmp(n, n);
  double dfv=0, out = 0;
  
#ifdef _OPENMP
  omp_set_dynamic(0);
#endif
  
//  double max_calc_value = -DBL_MAX; // minimum double value
#pragma omp parallel num_threads(nthreads)
{
  double dfv_thread, out_thread=0;
  
#pragma omp for
  for (int k = 1; k <= n; k++)
  {
//  R_CheckUserInterrupt();
    double dfvtmp;
    dfv_thread = 0;
    for (int i = 0; i < n; i++)
    {
      for (int j =0; j < n; j++)
      {
        if ((yrank[i] >= yrank[k - 1]) && (xrank[j] >= xrank[k - 1]))
          tmp(i, j) += 1;
        dfvtmp = tmp(i, j)/n - M(i, j)/n * k/n;
        if (fabs(dfvtmp) > dfv)
          dfv_thread = dfvtmp;
      }
    }
    if (dfv_thread > out_thread)
      out_thread = dfv_thread;
  }
  
#pragma omp critical
{
  if (out_thread > out) 
  {
    out = out_thread;
  }
}
}
  return out;
}



// ******************************************************************
