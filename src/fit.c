#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#include "dist.h"
#include "lambda.h"

// *******************************************************************************

double cloglkhd(double *tht, SEXP rdata, int verbose);
void cloglkhdGr(double *tht, SEXP rdata, int verbose, double *fv, double *dfv);

// *******************************************************************************

void clinesearch(SEXP rdata,
	    double *xOld,
	    double *h,
	    double *fv,
	    int verbose,
	    double *ram)
{
  R_CheckUserInterrupt();
  double const2 = 1.0e-16, ram1, ram2, ram3, fv1, fv2, fv3,
    xNew[8], a1, a2, a3, b1, b2;

  if (*ram <= 1.0e-30)
    *ram = 0.1;

  double hnorm = norm(h, 8);
  if (hnorm > 1)
    *ram = *ram/hnorm;

  ram1 = 0;
  ram2 = *ram;
  fv1  = *fv;

  for (int i = 0; i < 8; i++)
    xNew[i] = xOld[i] + ram2 * h[i];
  fv2 = cloglkhd(xNew, rdata, verbose);

  if (fv2 > fv1)
    goto stat50;

 stat30:
  ram3 = ram2*2.0;
  for (int i = 0; i < 8 ; i++)
    xNew[i] = xOld[i] + ram3 * h[i];
  fv3 = cloglkhd(xNew, rdata, verbose);
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
  for (int i = 0; i < 8; i++)
    xNew[i] = xOld[i] + ram2 * h[i];
  fv2 = cloglkhd(xNew, rdata, verbose);
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
      for (int i = 0; i < 8; i++)
	xNew[i] = xOld[i] + *ram*h[i];
      *fv = cloglkhd(xNew, rdata, verbose);
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
      for (int i = 0; i < 8; i++)
	xNew[i] = xOld[i] + *ram*h[i];
      *fv = cloglkhd(xNew, rdata, verbose);
      if (fv2 < *fv)
	*ram = ram2;
      return;
    }
}

// *******************************************************************************
// log-likelihood function of the model

double cloglkhd(double *tht,
                SEXP rdata,
                int verbose)
{
  // extract events
  SEXP revents = VECTOR_ELT(rdata, 0), dim;
  PROTECT(dim = allocVector(INTSXP, 2));
  dim = getAttrib(revents, R_DimSymbol);
  int N = INTEGER(dim)[0];
  double *events = REAL(revents);
  double t[N], x[N], y[N], m[N], bk[N];
  int flag[N];
  for (int i = 0; i < N; i++)
    {
      t[i] = events[i];
      x[i] = events[N + i];
      y[i] = events[2 * N + i];
      m[i] = events[3 * N + i];
      flag[i] = events[4 * N + i];
      bk[i] = events[5 * N + i];
    }

  // extract polygon information
  SEXP rpoly = VECTOR_ELT(rdata, 1), pdim;
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
  SEXP tperiod = VECTOR_ELT(rdata, 2);
  double *tper = REAL(tperiod);
  double tstart2 = tper[0], tlength = tper[1];

  // extract integral of spatial intensity over the obsevation window
  SEXP rinteg0 = VECTOR_ELT(rdata, 3);
  double integ0 = REAL(rinteg0)[0];

  double s , fv1 = 0, fv2 = 0, fv;

  for (int j = 0; j < N; j++)
    {
      if (flag[j] == 1)
	{
          s = clambdaj(tht, j, t, x, y, m, bk);
	  if (s > 1.0e-25)
	    fv1 += log(s);
	  else
	    fv1 -= 100.0;
	}
      fv2 += cintegj(tht, j, t, x, y, m, &np, px, py, &tstart2, &tlength);
    }

  fv2 += tht[0] * tht[0] * (integ0);

  fv = -fv1 + fv2;

  if (verbose == 1)
    Rprintf("Function Value = %8.5f\t%7.2f\t%7.2f\n", fv, -fv1, fv2);
  UNPROTECT(2);
  return fv;
}

// *******************************************************************************

void cloglkhdGr(double *tht,
	       SEXP rdata,
	       int verbose,
	       double *fv,
	       double *dfv)
{
  // extract events
  SEXP revents = VECTOR_ELT(rdata, 0), dim;
  PROTECT(dim = allocVector(INTSXP, 2));
  dim = getAttrib(revents, R_DimSymbol);
  int N = INTEGER(dim)[0];
  double *events = REAL(revents);
  double t[N], x[N], y[N], m[N], bk[N];
  int flag[N];
  for (int i = 0; i < N; i++)
    {
      t[i] = events[i];
      x[i] = events[N + i];
      y[i] = events[2 * N + i];
      m[i] = events[3 * N + i];
      flag[i] = events[4 * N + i];
      bk[i] = events[5 * N + i];
    }

  // extract polygon information
  SEXP rpoly = VECTOR_ELT(rdata, 1), pdim;
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
  SEXP tperiod = VECTOR_ELT(rdata, 2);
  double *tper = REAL(tperiod);
  double tstart2 = tper[0], tlength = tper[1];

  // extract integral of spatial intensity over the obsevation window
  SEXP rinteg0 = VECTOR_ELT(rdata, 3);
  double integ0 = REAL(rinteg0)[0];

  double fv1temp, g1temp[8], fv2temp, g2temp[8], fv1 = 0, fv2 = 0, df1[8] = {0}, df2[8] = {0};

  for (int j = 0; j < N; j++)
    {
      if (flag[j] == 1)
	{
          clambdajGr(tht, j, t, x, y, m, bk, &fv1temp, &g1temp[0]);

	  if (fv1temp > 1.0e-25)
	    fv1 += log(fv1temp);
	  else
	    fv1 -= 100.0;

	  for (int i = 0; i < 8; i++)
	    df1[i] += g1temp[i] / fv1temp;
	}

      cintegjGr(tht, j, t, x, y, m, &np, px, py, &tstart2, &tlength, &fv2temp, &g2temp[0]);
      fv2 += fv2temp;
      for (int i = 0; i < 8; i++)
	df2[i] += g2temp[i];
    }

  fv2 += tht[0] * tht[0] * (integ0);
  df2[0] = integ0 * tht[0] * 2;

  *fv = -fv1 + fv2;

  for (int i = 0; i < 8; i++)
    dfv[i] = -df1[i] + df2[i];
  if (verbose == 1)
    {
      Rprintf("Function Value = %8.2f\t%7.2f\t%7.2f\n", *fv, -fv1, fv2);
      for (int i = 0; i < 8; i++)
	Rprintf("Gradiant [%d] = %8.2f\ttheta[%d] = %2.8f\n", i + 1, dfv[i], i + 1, tht[i]);
    }
  UNPROTECT(2);
  return;
}

// *******************************************************************************

SEXP cfit(SEXP theta, SEXP rdata, SEXP ihess, SEXP rverbose)
{
  SEXP estimate, fvout, dfvout, aic, hess, out;

  // extract model parameters
  double *tht = REAL(theta);

  // extract the inverse of hessian matrix
  double *cihess = REAL(ihess);

  // extract verbose control
  int verbose = INTEGER(rverbose)[0];

  if (verbose == 1)
    Rprintf("\tstart Davidon-Fletcher-Powell procedure ... \n");

  double tau1 = 1.0e-6, tau2 = 1.0e-6, eps1 = 1.0e-6, eps2 = 1.0e-6,
    const1 = 1.0e-17;

  double ramda = 0.05, fv, sum, stem, s1, s2, ss;
  double h[8][8] = {{0}}, s[8] = {0}, dx[8] = {0}, g0[8] = {0}, g[8] = {0}, dg[8], wrk[8];

  // Initial estimate of inverse of hessian matrix
  for (int i = 0; i < 8; i++)
    for (int j = 0; j < 8; j++)
      h[i][j] = cihess[i * 8 + j];

  cloglkhdGr(tht, rdata, verbose, &fv, g);

  for (int iter = 1; iter < 10; iter++)
  {
    R_CheckUserInterrupt();
    for (int ic = 0; ic < 8; ic++)
    {
      if (ic > 0 || iter > 1)
      {
        for (int i = 0; i < 8; i++)
          dg[i] = g[i] - g0[i];

        for (int i = 0; i < 8; i++)
        {
          sum = 0.0;
          for (int j = 0; j < 8; j++)
            sum += dg[j] * h[i][j];
          wrk[i] = sum;
        }

        s1 = 0.0;
        s2 = 0.0;
        for (int i = 0; i < 8; i++)
        {
          s1 += wrk[i] * dg[i];
          s2 += dx[i] * dg[i];
        }

        if (s1 <= const1 || s2 <= const1)
        {
          PROTECT(estimate = allocVector(REALSXP, 8));
          PROTECT(fvout = allocVector(REALSXP, 1));
          PROTECT(dfvout = allocVector(REALSXP, 8));
          PROTECT(aic = allocVector(REALSXP, 1));
          PROTECT(hess = allocVector(REALSXP, 64));
          double *estimP = REAL(estimate), *fvP = REAL(fvout),
            *dfvP = REAL(dfvout), *aicP=REAL(aic), *hessP = REAL(hess);
          fvP[0] = -fv;
          aicP[0] = 2 * (fv + 8);
          if (verbose == 1)
            Rprintf ("loglikelihood = %8.5f\tAIC = %8.5f\n",
                     -fv, 2 * (fv + 8));
          for( int i = 0; i < 8; i++ )
          {
            dfvP[i] = g[i];
            estimP[i] = tht[i];
            for (int j = 0; j < 8; j++)
              hessP[i + 8 * j] = h[i][j];
            if (verbose == 1)
              Rprintf("theta[%d] = %2.8f\t gradient[%d] = %8.4f\n",
                      i + 1, pow(tht[i], 2), i + 1, g[i]);
          }
          PROTECT(out = allocVector(VECSXP, 5));
          SET_VECTOR_ELT(out, 0, estimate);
          SET_VECTOR_ELT(out, 1, fvout);
          SET_VECTOR_ELT(out, 2, dfvout);
          SET_VECTOR_ELT(out, 3, aic);
          SET_VECTOR_ELT(out, 4, hess);
          UNPROTECT(6);
          return(out);
        }

        if (s1 <= s2)
        {
          // fletcher type correction
          stem = s1 / s2 + 1.0;
          for (int i = 0; i < 8; i++)
            for (int j = i; j < 8; j++)
            {
              h[i][j] -= (dx[i] * wrk[j] + wrk[i] * dx[j] - dx[i] * dx[j] * stem) / s2;
              h[j][i] = h[i][j];
            }
        }
        else
        {
          // Update the inverse of hessian matrix
          for (int i = 0; i < 8; i++)
            for (int j = i; j < 8; j++)
            {	// davidon-fletcher-powell type correction
              h[i][j] += dx[i] * dx[j]/s2 - wrk[i] * wrk[j] / s1;
              h[j][i] = h[i][j];
            }
        }
      }
      ss = 0.0;
      for (int i = 0; i < 8; i++)
      {
        sum = 0.0;
        for (int j = 0; j < 8; j++)
          sum += h[i][j] * g[j];
        ss += sum * sum;
        s[i] = -sum;
      }
      s1 = 0.0;
      s2 = 0.0;
      for (int i = 0; i < 8; i++)
      {
        s1 += s[i] * g[i];
        s2 += g[i] * g[i];
      }

      double ds2 = sqrt(s2);
      double gtem = fabs(s1) / ds2;
      if (gtem <= tau1  &&  ds2 <= tau2)
      {
        PROTECT(estimate = allocVector(REALSXP, 8));
        PROTECT(fvout = allocVector(REALSXP, 1));
        PROTECT(dfvout = allocVector(REALSXP, 8));
        PROTECT(aic = allocVector(REALSXP, 1));
        PROTECT(hess = allocVector(REALSXP, 64));
        double *estimP = REAL(estimate), *fvP = REAL(fvout),
          *dfvP = REAL(dfvout), *aicP=REAL(aic), *hessP=REAL(hess);
        fvP[0] = -fv;
        aicP[0] = 2 * (fv + 8);
        if (verbose == 1)
          Rprintf ("loglikelihood = %8.5f\tAIC = %8.5f\n", -fv, 2 * (fv + 8));
        for( int i = 0; i < 8; i++ )
        {
          dfvP[i] = g[i];
          estimP[i] = tht[i];
          for (int j = 0; j < 8; j++)
            hessP[i + 8 * j] = h[i][j];
          if (verbose == 1)
            Rprintf("theta[%d] = %2.8f\t gradient[%d] = %8.4f\n", i + 1, pow(tht[i], 2), i + 1, g[i]);
        }
        PROTECT(out = allocVector(VECSXP, 5));
        SET_VECTOR_ELT(out, 0, estimate);
        SET_VECTOR_ELT(out, 1, fvout);
        SET_VECTOR_ELT(out, 2, dfvout);
        SET_VECTOR_ELT(out, 3, aic);
        SET_VECTOR_ELT(out, 4, hess);
        UNPROTECT(6);
        return(out);
      }

      if (s1 >= 0)
      {
        for (int i = 0; i < 8; i++)
        {
          for (int j = 0; j < 8; j++)
            h[i][j] = 0.0;
          h[i][i] = 1.0;
          s[i] = -s[i];
        }
      }


      double ed = fv;
      // line  search
      if (verbose == 1)
        Rprintf("\nstart line search along the specified direction ...\n");
      clinesearch(rdata, tht, s, &ed, verbose, &ramda);
      if (verbose == 1)
        Rprintf("back to Davidon-Fletcher-Powell Procedure: zeta = %f\n", ramda);

      R_CheckUserInterrupt();

      s1 = 0.0;
      for (int i = 0; i < 8; i++)
      {
        dx[i] = s[i] * ramda;
        s1 += dx[i] * dx[i];
        g0[i] = g[i];
        tht[i] += dx[i];
      }
      double fv0 = fv;
      cloglkhdGr(tht, rdata, verbose, &fv, g);

      s2 = 0.0;
      for (int i = 0; i < 8; i++)
        s2 += g[i] * g[i];
      if (sqrt(s2) > tau2)
        continue;
      if (fv0/fv - 1.0 < eps1 && sqrt(s1) < eps2)
      {
        PROTECT(estimate = allocVector(REALSXP, 8));
        PROTECT(fvout = allocVector(REALSXP, 1));
        PROTECT(dfvout = allocVector(REALSXP, 8));
        PROTECT(aic = allocVector(REALSXP, 1));
        PROTECT(hess = allocVector(REALSXP, 64));
        double *estimP = REAL(estimate), *fvP = REAL(fvout),
          *dfvP = REAL(dfvout), *aicP=REAL(aic), *hessP=REAL(hess);
        fvP[0] = -fv;
        aicP[0] = 2 * (fv + 8);
        if (verbose == 1)
          Rprintf ("loglikelihood = %8.5f\tAIC = %8.5f\n", -fv, 2 * (fv + 8));
        for( int i = 0; i < 8; i++ )
        {
          dfvP[i] = g[i];
          estimP[i] = tht[i];
          for (int j = 0; j < 8; j++)
            hessP[i + 8 * j] = h[i][j];
          if (verbose == 1)
            Rprintf("theta[%d] = %2.8f\t gradient[%d] = %8.4f\n", i + 1, pow(tht[i], 2), i + 1, g[i]);
        }
        PROTECT(out = allocVector(VECSXP, 5));
        SET_VECTOR_ELT(out, 0, estimate);
        SET_VECTOR_ELT(out, 1, fvout);
        SET_VECTOR_ELT(out, 2, dfvout);
        SET_VECTOR_ELT(out, 3, aic);
        SET_VECTOR_ELT(out, 4, hess);
        UNPROTECT(6);
        return(out);
      }
    }
  }
  return(0);
}

// *******************************************************************************

