#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "dist.h"
#include "poly.h"

// *******************************************************************************

double fr(double r, double w[])
{
  double gamma = w[0], D = w[1], q = w[2], mag = w[3], sig = D * exp(gamma * mag);
  return (1 - pow(1 + r * r / sig, 1 - q)) / (2 * PI);
}

double dgamma_fr(double r, double w[])
{
  double gamma = w[0], D = w[1], q = w[2], mag = w[3], sig = D * exp(gamma * mag);
  return (1 -q) * pow(1 + r * r / sig, -q) * mag * r * r / sig / (2 * PI);
}

double dD_fr(double r, double w[])
{
  double gamma = w[0], D = w[1], q = w[2], mag = w[3], sig = D * exp(gamma * mag);
  return (1 - q) * pow(1 + r * r / sig, -q) / D * r * r / sig / (2 * PI);
}

double dq_fr(double r, double w[])
{
  double gamma = w[0], D = w[1], q = w[2], mag = w[3], sig = D*exp( gamma*mag );
  return pow(1 + r * r / sig, 1 - q) * log(1 + r * r / sig) / (2 * PI);
}

// *******************************************************************************

double clambdaj(double *theta,
               int j,
               double *t,
               double *x,
               double *y,
               double *m,
               double *bk)
{
  // extract model parameters
  double mu = theta[0] * theta[0],
    A     = theta[1] * theta[1],
    c     = theta[2] * theta[2],
    alpha = theta[3] * theta[3],
    p     = theta[4] * theta[4],
    D     = theta[5] * theta[5],
    q     = theta[6] * theta[6],
    gamma = theta[7] * theta[7];

  double part1, part2, part3, delta, sig, r2;
  double s = mu * bk[j];

  for (register int i = 0; i < j; i++)
    {
      part1 = exp(alpha * m[i]);

      delta = t[j] - t[i];
      part2 = (p - 1)/c * pow(1 + delta/c, - p);

      sig = D * exp(gamma * m[i]);
      r2 = dist2(x[j], y[j], x[i], y[i]);
      part3 = (q - 1) / (sig * PI) * pow(1 + r2 / sig, - q);

      s += A * part1 * part2 * part3;
    }
  return s;
}

// *******************************************************************************

void clambdajGr(double *theta,
               int j,
               double *t,
               double *x,
               double *y,
               double *m,
               double *bk,
	       double *fv,
	       double *dfv)
{
  // extract model parameters
  double mu = theta[0] * theta[0],
    A     = theta[1] * theta[1],
    c     = theta[2] * theta[2],
    alpha = theta[3] * theta[3],
    p     = theta[4] * theta[4],
    D     = theta[5] * theta[5],
    q     = theta[6] * theta[6],
    gamma = theta[7] * theta[7];

  double part1, part2, part3, part1_alpha, part2_c, part2_p, part3_d, part3_q,
    part3_gamma, delta, sig, r2, sg1, sg2 = 0, sg3 = 0,
    sg4 = 0, sg5 = 0, sg6 = 0, sg7 = 0, sg8 = 0;

  double s = mu * bk[j];
  sg1 = bk[j];

  for (register int i = 0; i < j; i++)
    {
      part1 = exp(alpha * m[i]);

      delta = t[j] - t[i];
      part2 = (p - 1)/c * pow(1 + delta / c, - p);

      sig   = D * exp(gamma * m[i]);
      r2 = dist2(x[j], y[j], x[i], y[i]);
      part3 = (q - 1)/(sig * PI) * pow(1 + r2/sig, - q);

      s    += A * part1 * part2 * part3;
      sg2  += part1 * part2 * part3;

      part2_c = part2 * (-1/c - p/(c + delta) + p/c);
      sg3    += A * part1 * part2_c * part3;

      part1_alpha = part1 * m[i];
      sg4        += A * part1_alpha * part2 * part3;

      part2_p = part2 * (1/(p - 1) - log(1 + delta/c));
      sg5    += A * part1 * part2_p * part3;

      part3_d = part3 / D * (-1 + q * (1 - 1/(1 + r2/sig)));
      sg6    += A * part1 * part2 * part3_d;

      part3_q = part3 * (1/(q - 1) - log(1 + r2/sig));
      sg7    += A * part1 * part2 * part3_q;

      part3_gamma = part3 * (-m[i] + q * m[i] * (1 - 1/(1 + r2/sig)));
      sg8        += A * part1 * part2 * part3_gamma;
    }

  *fv      = s;
  dfv[ 0 ] = sg1 * 2 * theta[0];
  dfv[ 1 ] = sg2 * 2 * theta[1];
  dfv[ 2 ] = sg3 * 2 * theta[2];
  dfv[ 3 ] = sg4 * 2 * theta[3];
  dfv[ 4 ] = sg5 * 2 * theta[4];
  dfv[ 5 ] = sg6 * 2 * theta[5];
  dfv[ 6 ] = sg7 * 2 * theta[6];
  dfv[ 7 ] = sg8 * 2 * theta[7];
}

// *******************************************************************************

double cintegj(double *theta,
              int j,
              double *t,
              double *x,
              double *y,
              double *m,
              int *np,
              double *px,
              double *py,
              double *tstart2,
              double *tlength)
{

  // extract model parameters
  double //mu = theta[0] * theta[0],
    A     = theta[1] * theta[1],
    c     = theta[2] * theta[2],
    alpha = theta[3] * theta[3],
    p     = theta[4] * theta[4],
    D     = theta[5] * theta[5],
    q     = theta[6] * theta[6],
    gamma = theta[7] * theta[7];

  double ttemp, ttemp1, ttemp2, gi, gi1, gi2, w[4], si, sk;

  if (t[j] > *tstart2)
    {
      ttemp = *tlength - t[j];
      gi  = 1 - pow(1 + ttemp/c, 1 - p);
    }
  else
    {
      ttemp1 = *tstart2 - t[j];
      ttemp2 = *tlength - t[j];

      gi1  = 1 - pow(1 + ttemp1/c, 1 - p);
      gi2  = 1 - pow(1 + ttemp2/c, 1 - p);
      gi   = gi2 - gi1;
    }

  w[ 0 ] = gamma;
  w[ 1 ] = D;
  w[ 2 ] = q;
  w[ 3 ] = m[j];

  si = polyinteg(fr, w, np, px, py, x[j], y[j]);
  sk = A * exp(alpha * m[j]);
  return sk * gi * si;
}

// *******************************************************************************

void cintegjGr(double *theta,
              int j,
              double *t,
              double *x,
              double *y,
              double *m,
              int *np,
              double *px,
              double *py,
              double *tstart2,
              double *tlength,
              double *fv,
	      double *dfv)
{
  // extract model parameters
  double //mu = theta[0] * theta[0],
    A     = theta[1] * theta[1],
    c     = theta[2] * theta[2],
    alpha = theta[3] * theta[3],
    p     = theta[4] * theta[4],
    D     = theta[5] * theta[5],
    q     = theta[6] * theta[6],
    gamma = theta[7] * theta[7];

  double ttemp, ttemp1, ttemp2, gi, gi1, gi2, gic, gic1, gic2, gip, gip1,
    gip2, w[4], si, sid, siq, sigamma, sk;

  if (t[j] > *tstart2)
    {
      ttemp = *tlength - t[j];

      gi  = 1 - pow(1 + ttemp/c, 1 - p);
      gic = - (1 - gi) * (1 - p) * ( 1/(c + ttemp) - 1/c);
      gip = - (1 - gi) * (log(c) - log(c + ttemp));
    }
  else
    {
      ttemp1 = *tstart2 - t[j];
      ttemp2 = *tlength - t[j];

      gi1  = 1 - pow(1 + ttemp1/c, 1 - p);
      gi2  = 1 - pow(1 + ttemp2/c, 1 - p);
      gic1 = - (1 - gi1) * (1 - p) * (1/(c + ttemp1) - 1/c);
      gic2 = - (1 - gi2) * (1 - p) * (1/(c + ttemp2) - 1/c);
      gip1 = - (1 - gi1) * (log(c) - log(c + ttemp1));
      gip2 = - (1 - gi2) * (log(c) - log(c + ttemp2));

      gi  = gi2 - gi1;
      gic = gic2 - gic1;
      gip = gip2 - gip1;
    }

  w[0] = gamma;
  w[1] = D;
  w[2] = q;
  w[3] = m[j];

  si      = polyinteg(fr, w, np, px, py, x[j], y[j]);
  sid     = polyinteg(dD_fr, w, np, px, py, x[j], y[j]);
  siq     = polyinteg(dq_fr, w, np, px, py, x[j], y[j]);
  sigamma = polyinteg(dgamma_fr, w, np, px, py, x[j], y[j]);

  sk = A * exp(alpha * m[j]);
  *fv      = sk * gi * si;
  dfv[ 0 ] = 0;
  dfv[ 1 ] = sk * gi  * si / A        * 2 * theta[1];
  dfv[ 2 ] = sk * gic * si            * 2 * theta[2];
  dfv[ 3 ] = sk * gi  * si * m[j]     * 2 * theta[3];
  dfv[ 4 ] = sk * gip * si            * 2 * theta[4];
  dfv[ 5 ] = sk * gi  * sid           * 2 * theta[5];
  dfv[ 6 ] = sk * gi  * siq           * 2 * theta[6];
  dfv[ 7 ] = sk * gi  * sigamma       * 2 * theta[7];
  return;
}

// *******************************************************************************

SEXP clambdax(SEXP rt,
             SEXP rx,
             SEXP ry,
             SEXP theta,
             SEXP revents)
{
  SEXP dim;

  // extract events
  PROTECT(dim = allocVector(INTSXP, 2));
  dim = getAttrib(revents, R_DimSymbol);
  int N = INTEGER(dim)[0];
  double *events = REAL(revents);
  double t[N], x[N], y[N], m[N];
  for (int i = 0; i < N; i++)
    {
      t[i] = events[i];
      x[i] = events[N + i];
      y[i] = events[2 * N + i];
      m[i] = events[3 * N + i];
    }

  // extract model parameters
  double *tht = REAL(theta);
  double //mu = tht[0] * tht[0],
    A     = tht[1] * tht[1],
    c     = tht[2] * tht[2],
    alpha = tht[3] * tht[3],
    p     = tht[4] * tht[4],
    D     = tht[5] * tht[5],
    q     = tht[6] * tht[6],
    gamma = tht[7] * tht[7];

  // extract arguments
  double tt = *REAL(rt), xx = *REAL(rx), yy = *REAL(ry);

  double part1, part2, part3, delta, sig, r2;

  double s = 0;

  int i = 0;
  while (t[i] < tt && i < N)
    {
      part1 = exp(alpha * m[i]);

      delta = tt - t[i];
      part2 = (p - 1)/c * pow(1 + delta/c, - p);

      sig = D * exp(gamma * m[i]);
      r2 = dist2(xx, yy, x[i], y[i]);
      part3 = (q - 1) / (sig * PI) * pow(1 + r2 / sig, - q);

      s += A * part1 * part2 * part3;
      i++;
    }

  SEXP out;
  PROTECT(out = allocVector(REALSXP, 1));
  double *outP = REAL(out);
  *outP = s;
  UNPROTECT(2);
  return out;
}

// *******************************************************************************

