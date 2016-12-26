#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "dist.h"
#include "poly.h"
#include "lambda.h"

// *******************************************************************************

// value of the gaussian kernel function at r with bandwidth sig
double dGauss(double r, double sig)
{
  return exp(-(r * r) /(2 * sig * sig)) / (2 * PI * sig * sig);
}

// integral of the gaussian kernel with bandwidth w[0] from 0 to r
double pGauss(double r, double w[])
{
  return (1 - exp(-(r * r) / (2 * w[0] * w[0]))) / (2 * PI);
}

// *******************************************************************************

SEXP cdeclust(SEXP theta,
             SEXP rbwd,
             SEXP revents,
             SEXP rpoly,
             SEXP tperiod)
{
  SEXP dim, pdim, out, integ0;

  // extract events
  PROTECT(dim = allocVector(INTSXP, 2));
  dim = getAttrib(revents, R_DimSymbol);
  int N = INTEGER(dim)[0];
  double *events = REAL(revents);
  double t[N], x[N], y[N], m[N], bk[N], pb[N], lam[N];
  for (int i = 0; i < N; i++)
    {
      t[i] = events[i];
      x[i] = events[N + i];
      y[i] = events[2 * N + i];
      m[i] = events[3 * N + i];
      bk[i] = events[5 * N + i];
      pb[i] = events[6 * N + i];
      lam[i] = events[7 * N + i];
    }

  // extract polygon information
  PROTECT(pdim = allocVector(INTSXP, 2));
  pdim = getAttrib(rpoly, R_DimSymbol);
  int np = INTEGER(pdim)[0];
  double *poly = REAL(rpoly);
  double px[np], py[np];
  for (int i = 0; i < np; i++)
    {
      px[i] = poly[i];
      py[i] = poly[np + i];
    }

  // extract time period information
  double *tper = REAL(tperiod);
  double tstart2 = tper[0], tlength = tper[1];

  // extract bandwidthes
  double *bwd = REAL(rbwd);

  // extract model paramters
  double *tht = REAL(theta);

  double s, r0, w[1];

  for (int i = 0; i < N; i++)
    {
      s = 0;
      for (int j = 0; j < N; j++)
	{
	  r0 = dist(x[i], y[i], x[j], y[j]);
	  s += pb[j] * dGauss(r0, bwd[j]);
	}
      bk[i] = s / (tlength - tstart2);
      events[5 * N + i] = bk[i];
    }

  s = 0;
  for (int i = 0; i < N; i++)
    {
      w[0] = bwd[i];
      s += pb[i] * polyinteg(pGauss, w, &np, px, py, x[i], y[i]);
      lam[i] = clambdaj(tht,i, t, x, y, m, bk);
      events[6 * N + i] = (tht[0] * tht[0] * bk[i]) / lam[i];
      events[7 * N + i] = lam[i];
    }

  PROTECT(out = allocVector(VECSXP, 2));
  PROTECT(integ0 = allocVector(REALSXP, 1));
  double *integ0P = REAL(integ0);
  integ0P[0] = s;
  SET_VECTOR_ELT(out, 0, revents);
  SET_VECTOR_ELT(out, 1, integ0);
  UNPROTECT(4);
  return(out);
}

// *******************************************************************************

/*
SEXP bkgd2(SEXP rpb,
           SEXP rbwd,
           SEXP rdata)
{

  // extract data
  SEXP revents = VECTOR_ELT(rdata, 0);
  SEXP rpoly = VECTOR_ELT(rdata, 1);
  SEXP tperiod = VECTOR_ELT(rdata, 4);
  SEXP dim, pdim;

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

  // extract polygon information
  PROTECT(pdim = allocVector(INTSXP, 2));
  pdim = getAttrib(rpoly, R_DimSymbol);
  int np = INTEGER(pdim)[0];
  double *poly = REAL(rpoly);
  double px[np], py[np];
  for (int i = 0; i < np; i++)
    {
      px[i] = poly[i];
      py[i] = poly[np + i];
    }

  // extract time period information
  double *tper = REAL(tperiod);
  double tstart2 = tper[0], tlength = tper[1];

  // extract bandwidthes
  double *bwd = REAL(rbwd);

  // extract bavkground probabilities
  double *pb = REAL(rpb);

  SEXP bk;
  PROTECT(bk = allocVector(REALSXP, N));
  double *bkP = REAL(bk);
  double s, r0, w[1], lam[N];

  for (int i = 0; i < N; i++)
    {
      s = 0;
      for (int j = 0; j < N; j++)
	{
	  r0 = dist(x[i], y[i], x[j], y[j]);
	  s += pb[j] * dGauss(r0, bwd[j]);
	}
      bkP[i] = s / (tlength - tstart2);
    }

  s = 0;
  for (int i = 0; i < N; i++)
    {
      w[0] = bwd[i];
      s += pb[i] * polyinteg(pGauss, w, &np, px, py, x[i], y[i]);
    }

  SEXP out, integ0;
  PROTECT(out = allocVector(VECSXP, 2));
  PROTECT(integ0 = allocVector(REALSXP, 1));
  double *integ0P = REAL(integ0);
  *integ0P = s;
  SET_VECTOR_ELT(out, 0, bk);
  SET_VECTOR_ELT(out, 1, integ0);
  UNPROTECT(5);
  return(out);
}

// *******************************************************************************


SEXP probs(SEXP theta,
           SEXP rbk,
           SEXP rbwd,
           SEXP rpb,
           SEXP rdata)
{
  // extract model parameters
  double *tht = REAL(theta);
  double mu = tht[0] * tht[0],
    A     = tht[1] * tht[1],
    c     = tht[2] * tht[2],
    alpha = tht[3] * tht[3],
    p     = tht[4] * tht[4],
    D     = tht[5] * tht[5],
    q     = tht[6] * tht[6],
    gamma = tht[7] * tht[7];

  // extract data
  SEXP revents = VECTOR_ELT(rdata, 0);
  SEXP rpoly = VECTOR_ELT(rdata, 1);
  SEXP tperiod = VECTOR_ELT(rdata, 4);
  SEXP dim, pdim;

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

  // extract polygon information
  PROTECT(pdim = allocVector(INTSXP, 2));
  pdim = getAttrib(rpoly, R_DimSymbol);
  int np = INTEGER(pdim)[0];
  double *poly = REAL(rpoly);
  double px[np], py[np];
  for (int i = 0; i < np; i++)
    {
      px[i] = poly[i];
      py[i] = poly[np + i];
    }

  // extract time period information
  double *tper = REAL(tperiod);
  double tstart2 = tper[0], tlength = tper[1];

  // extract bandwidthes,
  double *bwd = REAL(rbwd);
  double *bk = REAL(rbk);
  double *pb = REAL(rpb);

  SEXP out, clam, integ0;
  PROTECT(out = allocVector(VECSXP, 3));
  PROTECT(clam = allocVector(REALSXP, N));
  PROTECT(integ0 = allocVector(REALSXP, 1));
  double *integ0P=REAL(integ0), *clamP=REAL(clam);

  double s=0, r0, w[1];

  s = 0;
  for (int i = 0; i < N; i++)
    {
      w[0] = bwd[i];
      s += pb[i] * polyinteg(pGauss, w, &np, px, py, x[i], y[i]);
      clamP[i] = lambdaj(tht,i, rdata);
      pb[i] = (mu * bk[i]) / clamP[i];
    }
//  SEXP cbk, cpb, clam, cinteg0;

  integ0P[0] = s;
  SET_VECTOR_ELT(out, 0, rpb);
  SET_VECTOR_ELT(out, 1, clam);
  SET_VECTOR_ELT(out, 2, integ0);
  UNPROTECT(5);
  return(out);
}

*/

