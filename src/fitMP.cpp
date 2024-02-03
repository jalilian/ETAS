#include <cmath>
#include <Rcpp.h>

#include "funcs.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// ******************************************************************
// approximating the integral of a function on a polygon region
// ******************************************************************

double polyintegXX(double (*func)(double, double, double []),
                   double m,
                   double funcpara[],
                                  NumericVector px,
                                  NumericVector py,
                                  double cx,
                                  double cy,
                                  int ndiv)
{
  int id;
  double sum = 0, dxx, dyy, x1, x2, y1, y2;
  double det, r0, r1, r2, theta;
  
  for (int k = 0; k < (px.length() - 1); ++k)
  {
    dxx = (px[k + 1] - px[k]) / ndiv;
    dyy = (py[k + 1] - py[k]) / ndiv;
    for (int l = 0; l < ndiv; ++l)
    {
      x1 = px[k] + dxx * l;
      y1 = py[k] + dyy * l;
      x2 = px[k] + dxx * (l + 1);
      y2 = py[k] + dyy * (l + 1);
      det = (x1 * y2 + y1 * cx + x2 * cy) - (x2 * y1 + y2 * cx + x1 * cy);
      
      if (fabs(det) < 1.0e-10)
        continue;
      
      id = 1;
      if (det < 0)
        id = -1;
      
      r1 = dist(x1, y1, cx, cy);
      r2 = dist(x2, y2, cx, cy);
      theta = (r1 * r1 + r2 * r2 - dist2(x1, y1, x2, y2))/(2 * r1 * r2);
      if (fabs(theta) > 1)
        theta = 1 - 1.0e-10;
      
      theta = acos(theta);
      
      if (r1 + r2 > 1.0e-20)
      {
        r0 = dist(x1 + r1/(r1 + r2) * (x2 - x1),
                  y1 + r1/(r1 + r2) * (y2 - y1), cx, cy);
        
        sum += id * (func(r1, m, funcpara) / 6 +
          func(r0, m, funcpara) * 2 / 3 +
          func(r2, m, funcpara) / 6) * theta;
      }
    }
  }
  
  return sum;
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
  
public:
  void set(NumericMatrix revents,
           NumericMatrix rpoly,
           NumericVector tperiod,
           double rinteg0,
           int rndiv);
  double mloglik(NumericVector theta);
  void mloglikGr(NumericVector theta,
                 double *fv,
                 double *df);
  void linesearch(NumericVector xOld,
                  double *h,
                  double *fv,
                  double *ram);
  List fitfun(NumericVector tht,
              NumericMatrix ihess,
              double eps,
              bool verbose);
  double mloglikMP(NumericVector theta,
                   int nthreads);
  void mloglikGrMP(NumericVector theta,
                   double *fv,
                   double *df,
                   int nthreads);
  void linesearchMP(NumericVector xOld,
                    double *h,
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
               int rndiv)
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
}

// ******************************************************************
// minus log likelihood function
// ******************************************************************

double etas::mloglik(NumericVector theta)
{
  const double mu = theta[0] * theta[0];
  const double A = theta[1] * theta[1];
  const double c = theta[2] * theta[2];
  const double alpha = theta[3] * theta[3];
  const double p = theta[4] * theta[4];
  const double D = theta[5] * theta[5];
  const double q= theta[6] * theta[6];
  const double gamma = theta[7] * theta[7];

  double kparam[] = {A, alpha};
  double gparam[] = {c, p};
  double fparam[] = {D, gamma, q};
  
  double fv1 = 0, fv2 = 0, si, gi;
  
  for (int j = 0; j < t.length(); ++j)
  {
    if (flag[j] == 1)
    {
      double sumpart = mu * bk[j];
      for (int i = 0; i < j; i++)
      {
        sumpart += kappafun(m[i], kparam) *
          gfun(t[j] - t[i], gparam) *
          ffun(dist2(x[j], y[j], x[i], y[i]), m[i], fparam);
      }
      
      if (sumpart > 1.0e-25)
        fv1 += log(sumpart);
      else
        fv1 += -100;
    }
    
    if (t[j] > tstart2)
    {
      gi = gfunint(tlength - t[j], gparam);
    }
    else
    {
      gi = gfunint(tlength - t[j], gparam) - gfunint(tstart2 - t[j], gparam);
    }
    
    si = 0;
    double dpx, dpy, x1, x2, y1, y2, det, r0, r1, r2, phi;
    for (int k = 0; k < (px.length() - 1); ++k)
    {
      dpx = (px[k + 1] - px[k]) / ndiv;
      dpy = (py[k + 1] - py[k]) / ndiv;
      for (int l = 0; l < ndiv; ++l)
      {
        x1 = px[k] + dpx * l;
        y1 = py[k] + dpy * l;
        x2 = px[k] + dpx * (l + 1);
        y2 = py[k] + dpy * (l + 1);
        
        det = (x1 * y2 + y1 * x[j] + x2 * y[j]) -
          (x2 * y1 + y2 * x[j] + x1 * y[j]);
        if (fabs(det) < 1.0e-10)
          continue;
        
        r1 = dist(x1, y1, x[j], y[j]);
        r2 = dist(x2, y2, x[j], y[j]);
        phi = (r1 * r1 + r2 * r2 - dist2(x1, y1, x2, y2))/(2 * r1 * r2);
        if (fabs(phi) > 1)
          phi = 1 - 1.0e-10;
        
        phi = acos(phi);
        
        if (r1 + r2 > 1.0e-20)
        {
          r0 = dist(x1 + r1/(r1 + r2) * (x2 - x1),
                    y1 + r1/(r1 + r2) * (y2 - y1), x[j], y[j]);
          
          si += sgn(det) * (frfunint(r1, m[j], fparam) / 6 +
            frfunint(r0, m[j], fparam) * 2 / 3 +
            frfunint(r2, m[j], fparam) / 6) * phi;
        }
      }
    }
    
    fv2 += kappafun(m[j], kparam) * gi * si;
  }
  
  fv2 += mu * integ0;
  
  return -fv1 + fv2;
}

// ******************************************************************
// gradient of minus log likelihood function
// ******************************************************************

void etas::mloglikGr(NumericVector theta,
                     double *fv,
                     double *dfv)
{
  const double mu = theta[0] * theta[0];
  const double A = theta[1] * theta[1];
  const double c = theta[2] * theta[2];
  const double alpha = theta[3] * theta[3];
  const double p = theta[4] * theta[4];
  const double D = theta[5] * theta[5];
  const double q= theta[6] * theta[6];
  const double gamma = theta[7] * theta[7];

  double kparam[] = {A, alpha};
  double gparam[] = {c, p};
  double fparam[] = {D, gamma, q};

  double fv1 = 0, fv2 = 0, df1[8] = {0}, df2[8] = {0};
  
  NumericVector part1(3), int_part1(3), part2(3), int_part2(3), part3(3);
  
  for (int j = 0; j < N; ++j)
  {
    if (flag[j] == 1)
    {
      double fv1temp = mu * bk[j];
      double g1temp[8] = {0};
      g1temp[0] = bk[j];
      
      for (int i = 0; i < j; i++)
      {
        part1 = dkappafun(m[i], kparam);
        
        part2 = dgfun(t[j] - t[i], gparam);
        
        part3 = dffun(dist2(x[j], y[j], x[i], y[i]), m[i], fparam);

        fv1temp    += part1[0] * part2[0] * part3[0];

        // part1_A
        g1temp[1]  += part1[1] * part2[0] * part3[0];
        
        // part2_csig
        g1temp[2] += part1[0] * part2[1] * part3[0];
        
        //part1_alpha
        g1temp[3]  += part1[2] * part2[0] * part3[0];
        
        // part2_p
        g1temp[4] += part1[0] * part2[2] * part3[0];
        
        // part3_d
        g1temp[5] += part1[0] * part2[0] * part3[1];
        
        // part3_q
        g1temp[6] += part1[0] * part2[0] * part3[2];
        
        // part3_gamma
        g1temp[7]  += part1[0] * part2[0] * part3[3];
      }
      
      if (fv1temp > 1.0e-25)
        fv1 += log(fv1temp);
      else
        fv1 += -100;
      
      for (int i = 0; i < 8; i++)
      {
        g1temp[i] *= 2 * theta[i];
        df1[i] += g1temp[i] / fv1temp;
      }
    }
    
    if (t[j] > tstart2)
    {
      int_part2 = dgfunint(tlength - t[j], gparam);
    }
    else
    {
      int_part2 = dgfunint(tlength - t[j], gparam) -
        dgfunint(tstart2 - t[j], gparam);
    }
    
    NumericVector int_part3(4);
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
          
          int_part3 += id * (dfrfunint(r1, m[j], fparam) / 6 +
              dfrfunint(r0, m[j], fparam) * 2.0 / 3 +
              dfrfunint(r2, m[j], fparam) / 6) * phi;
        }
      }
    }
    
    int_part1 = dkappafun(m[j], kparam);

    double fv2temp  = int_part1[0] * int_part2[0] * int_part3[0];
    double g2temp[8] = {0};

    // d A
    g2temp[ 1 ] = int_part1[1] * int_part2[0]  * int_part3[0];
    // d c
    g2temp[ 2 ] = int_part1[0] * int_part2[1] * int_part3[0];
    // d alpha
    g2temp[ 3 ] = int_part1[2] * int_part2[0]  * int_part3[0];
    // d p
    g2temp[ 4 ] = int_part1[0] * int_part2[2] * int_part3[0];
    // d D
    g2temp[ 5 ] = int_part1[0] * int_part2[0]  * int_part3[1];
    // d q
    g2temp[ 6 ] = int_part1[0] * int_part2[0]  * int_part3[2];
    // d gamma
    g2temp[ 7 ] = int_part1[0] * int_part2[0]  * int_part3[3];
    
    fv2 += fv2temp;
    for (int i = 0; i < 8; i++)
    {
      g2temp[i] *= 2 * theta[i];
      df2[i] += g2temp[i];
    }
  }
  
  fv2 += mu * integ0;
  df2[0] = integ0 * theta[0] * 2;
  
  *fv = -fv1 + fv2;
  
  for (int i = 0; i < 8; ++i)
    dfv[i] = -df1[i] + df2[i];
  return;
}


// ******************************************************************
// line search for the optimization algorithm
// ******************************************************************

void etas::linesearch(NumericVector xOld,
                      double *h,
                      double *fv,
                      double *ram)
{
  R_CheckUserInterrupt();
  double const2 = 1.0e-16, ram1, ram2, ram3, fv1, fv2, fv3,
    a1, a2, a3, b1, b2;
  
  NumericVector xNew(8);
  
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
  fv2 = mloglik(xNew);
  
  if (fv2 > fv1)
    goto stat50;
  
  stat30:
    ram3 = ram2*2.0;
  for (int i = 0; i < 8 ; i++)
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
  for (int i = 0; i < 8; i++)
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
    for (int i = 0; i < 8; i++)
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
    for (int i = 0; i < 8; i++)
      xNew[i] = xOld[i] + *ram*h[i];
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
  NumericVector estimate(8), dfvout(8);
  double fvout, aic;
  
  if (verbose)
    Rprintf("\tstart Davidon-Fletcher-Powell procedure ... \n");
  
  double tau1 = eps, tau2 = eps, eps1 = eps, eps2 = eps,
    const1 = 1.0e-17;
  
  double ramda = 0.05, fv, s1, s2;
  double h[8][8], s[8] = {0}, dx[8] = {0}, g0[8] = {0},
    g[8] = {0}, dg[8], wrk[8];
  
  // Initial estimate of inverse of hessian matrix
  for (int i = 0; i < 8; i++)
    for (int j = 0; j < 8; j++)
      h[i][j] = ihess(i, j);
  
  mloglikGr(tht, &fv, g);
  
  if (verbose)
  {
    Rprintf("Function Value = %8.4f\n", fv);
    for (int i = 0; i < 8; ++i)
      Rprintf("Gradient[%d] = %8.2f\ttheta[%d] = %2.6f\n", i + 1,
              g[i], i + 1, tht[i]);
  }
  
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
          double sum = 0;
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
          fvout = -fv;
          aic = 2 * (fv + 8);
          if (verbose)
            Rprintf ("loglikelihood = %8.5f\tAIC = %8.5f\n", -fv, 2 * (fv + 8));
          for( int i = 0; i < 8; i++ )
          {
            dfvout[i] = g[i];
            estimate[i] = tht[i];
            for (int j = 0; j < 8; j++)
              ihess(i, j) = h[i][j];
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
          for (int i = 0; i < 8; i++)
            for (int j = i; j < 8; j++)
            {
              h[i][j] -= (dx[i] * wrk[j] + wrk[i] * dx[j] -
                dx[i] * dx[j] * (1 + s1 / s2)) / s2;
              h[j][i] = h[i][j];
            }
        }
        else
        {
          // Update the inverse of hessian matrix
          for (int i = 0; i < 8; i++)
            for (int j = i; j < 8; j++)
            {	// davidon-fletcher-powell type correction
              h[i][j] += dx[i] * dx[j] / s2 - wrk[i] * wrk[j] / s1;
              h[j][i] = h[i][j];
            }
        }
      }
      
      double ss = 0;
      for (int i = 0; i < 8; i++)
      {
        double sum = 0;
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
      
      if ((fabs(s1) / sqrt(s2) <= tau1) && (sqrt(s2) <= tau2))
      {
        fvout = -fv;
        aic = 2 * (fv + 8);
        if (verbose)
          Rprintf ("loglikelihood = %8.5f\tAIC = %8.5f\n", -fv, 2 * (fv + 8));
        for( int i = 0; i < 8; i++ )
        {
          dfvout[i] = g[i];
          estimate[i] = tht[i];
          for (int j = 0; j < 8; j++)
            ihess(i, j) = h[i][j];
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
        for (int i = 0; i < 8; i++)
        {
          for (int j = 0; j < 8; j++)
            h[i][j] = 0.0;
          h[i][i] = 1.0;
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
      for (int i = 0; i < 8; i++)
      {
        dx[i] = s[i] * ramda;
        s1 += dx[i] * dx[i];
        g0[i] = g[i];
        tht[i] += dx[i];
      }
      
      double fv0 = fv;
      mloglikGr(tht, &fv, g);
      
      if (verbose)
      {
        Rprintf("Function Value = %8.4f\n", fv);
        for (int i = 0; i < 8; ++i)
          Rprintf("Gradient[%d] = %8.2f\ttheta[%d] = %2.6f\n", i + 1,
                  g[i], i + 1, tht[i]);
      }
      
      s2 = 0;
      for (int i = 0; i < 8; i++)
        s2 += g[i] * g[i];
      if (sqrt(s2) > tau2)
        continue;
      
      if (fv0/fv - 1 < eps1 && sqrt(s1) < eps2)
      {
        fvout = -fv;
        aic = 2 * (fv + 8);
        if (verbose)
          Rprintf ("loglikelihood = %8.5f\tAIC = %8.5f\n", -fv, 2 * (fv + 8));
        for( int i = 0; i < 8; i++ )
        {
          dfvout[i] = g[i];
          estimate[i] = tht[i];
          for (int j = 0; j < 8; j++)
            ihess(i, j) = h[i][j];
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
  const double mu = theta[0] * theta[0];
  const double A = theta[1] * theta[1];
  const double c = theta[2] * theta[2];
  const double alpha = theta[3] * theta[3];
  const double p = theta[4] * theta[4];
  const double D = theta[5] * theta[5];
  const double q= theta[6] * theta[6];
  const double gamma = theta[7] * theta[7];

  double kparam[] = {A, alpha};
  double gparam[] = {c, p};
  double fparam[] = {D, gamma, q};

  double fv1 = 0, fv2 = 0;
  
#pragma omp parallel num_threads(nthreads)
{
  double fv1_thread = 0, fv2_thread = 0;
  
#pragma omp for
  for (int j = 0; j < N; ++j)
  {
    double s_thread, gi;
    if (flag[j] == 1)
    {
      s_thread = mu * bk[j];
      for (int i = 0; i < j; ++i)
      {
        s_thread += kappafun(m[i], kparam) *
          gfun(t[j] - t[i], gparam) *
          ffun(dist2(x[j], y[j], x[i], y[i]), m[i], fparam);
      }
      
      if (s_thread > 1.0e-25)
        fv1_thread += log(s_thread);
      else
        fv1_thread += -100;
    }
    
    if (t[j] > tstart2)
    {
      gi = gfunint(tlength - t[j], gparam);
    }
    else
    {
      gi = gfunint(tlength - t[j], gparam) - gfunint(tstart2 - t[j], gparam);
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
                    y1 + r1/(r1 + r2) * (y2 - y1),
                    x[j], y[j]);
          
          si += sgn(det) * (frfunint(r1, m[j], fparam) / 6 +
            (frfunint(r0, m[j], fparam) * 2) / 3 +
            frfunint(r2, m[j], fparam) / 6) * phi;
        }
      }
    }
    
    fv2_thread += kappafun(m[j], kparam) * gi * si;
  }
#pragma omp critical
{
  fv1 += fv1_thread;
  fv2 += fv2_thread;
}
}

fv2 += mu * integ0;
return -fv1 + fv2;
}

// ******************************************************************
// gradient of minus log likelihood function: parallel computing
// ******************************************************************

void etas::mloglikGrMP(NumericVector theta,
                       double *fv,
                       double *dfv,
                       int nthreads)
{
  const double mu = theta[0] * theta[0];
  const double A = theta[1] * theta[1];
  const double c = theta[2] * theta[2];
  const double alpha = theta[3] * theta[3];
  const double p = theta[4] * theta[4];
  const double D = theta[5] * theta[5];
  const double q= theta[6] * theta[6];
  const double gamma = theta[7] * theta[7];

  double kparam[] = {A, alpha};
  double gparam[] = {c, p};
  double fparam[] = {D, gamma, q};

  double fv1 = 0, fv2 = 0, df1[8] = {0}, df2[8] = {0};
  
#pragma omp parallel num_threads(nthreads)
{
  double fv1_thread = 0, fv2_thread = 0,
    df1_thread[8] = {0}, df2_thread[8] = {0};
  
#pragma omp for //schedule(static)
  for (int j = 0; j < N; ++j)
  {
    if (flag[j] == 1)
    {
      double fv1temp = mu * bk[j];
      double g1temp[8] = {0};
      g1temp[0] = bk[j];
      
      g1temp[1] = g1temp[2] = g1temp[3] = g1temp[4] = 0;
      g1temp[5] = g1temp[6] = g1temp[7] = 0;
      
      for (int i = 0; i < j; i++)
      {
        NumericVector part1 = dkappafun(m[i], kparam);
        NumericVector part2 = dgfun(t[j] - t[i], gparam);
        NumericVector part3 = dffun(dist2(x[j], y[j], x[i], y[i]), m[i], fparam);

        fv1temp    += part1[0] * part2[0] * part3[0];

        // part1_A
        g1temp[1]  += part1[1] * part2[0] * part3[0];

        // part2_csig
        g1temp[2] += part1[0] * part2[1] * part3[0];

        //part1_alpha
        g1temp[3]  += part1[2] * part2[0] * part3[0];

        // part2_p
        g1temp[4] += part1[0] * part2[2] * part3[0];

        // part3_d
        g1temp[5] += part1[0] * part2[0] * part3[1];

        // part3_q
        g1temp[6] += part1[0] * part2[0] * part3[2];

        // part3_gamma
        g1temp[7]  += part1[0] * part2[0] * part3[3];
      }
      
      if (fv1temp > 1.0e-25)
        fv1_thread += log(fv1temp);
      else
        fv1_thread += -100;
      
      for (int i = 0; i < 8; i++)
      {
        g1temp[i] *= 2 * theta[i];
        df1_thread[i] += g1temp[i] / fv1temp;
      }
    }
    
    if (t[j] > tstart2)
    {
      NumericVector int_part2 = dgfunint(tlength - t[j], gparam);
    }
    else
    {
      NumericVector int_part2 = dgfunint(tlength - t[j], gparam) -
        dgfunint(tstart2 - t[j], gparam);
    }

    NumericVector int_part3(4);
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
                    y1 + r1/(r1 + r2) * (y2 - y1),
                    x[j], y[j]);
          
          int_part3 += id * (dfrfunint(r1, m[j], fparam) / 6 +
              dfrfunint(r0, m[j], fparam) * 2.0 / 3 +
              dfrfunint(r2, m[j], fparam) / 6) * phi;
        }
      }
    }

    NumericVector int_part1 = dkappafun(m[j], kparam);

    double fv2temp  = int_part1[0] * int_part2[0] * int_part3[0];
    double g2temp[8] = {0};

    // d A
    g2temp[ 1 ] = int_part1[1] * int_part2[0]  * int_part3[0];
    // d c
    g2temp[ 2 ] = int_part1[0] * int_part2[1] * int_part3[0];
    // d alpha
    g2temp[ 3 ] = int_part1[2] * int_part2[0]  * int_part3[0];
    // d p
    g2temp[ 4 ] = int_part1[0] * int_part2[2] * int_part3[0];
    // d D
    g2temp[ 5 ] = int_part1[0] * int_part2[0]  * int_part3[1];
    // d q
    g2temp[ 6 ] = int_part1[0] * int_part2[0]  * int_part3[2];
    // d gamma
    g2temp[ 7 ] = int_part1[0] * int_part2[0]  * int_part3[3];
    
    fv2_thread += fv2temp;
    for (int i = 0; i < 8; i++)
    {
      g2temp[i] *= 2 * theta[i];
      df2_thread[i] += g2temp[i];
    }
  }
#pragma omp critical
{
  fv1 += fv1_thread;
  fv2 += fv2_thread;
  for (int i = 0; i < 8; ++i)
  {
    df1[i] += df1_thread[i];
    df2[i] += df2_thread[i];
  }
}
}
fv2 += mu * integ0;
df2[0] = integ0 * theta[0] * 2;

*fv = -fv1 + fv2;
for (int i = 0; i < 8; ++i)
  dfv[i] = -df1[i] + df2[i];

return;
}


// ******************************************************************
// line search for the optimization algorithm: parallel computing
// ******************************************************************

void etas::linesearchMP(NumericVector xOld,
                        double *h,
                        double *fv,
                        double *ram,
                        int nthreads)
{
  R_CheckUserInterrupt();
  double const2 = 1.0e-16, ram1, ram2, ram3, fv1, fv2, fv3,
    a1, a2, a3, b1, b2;
  
  NumericVector xNew(8);
  
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
  fv2 = mloglikMP(xNew, nthreads);
  
  if (fv2 > fv1)
    goto stat50;
  
  stat30:
    ram3 = ram2*2.0;
  for (int i = 0; i < 8 ; i++)
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
  for (int i = 0; i < 8; i++)
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
    for (int i = 0; i < 8; i++)
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
    for (int i = 0; i < 8; i++)
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
  NumericVector estimate(8), dfvout(8);
  double fvout, aic;
  
  if (verbose)
    Rprintf("\tstart Davidon-Fletcher-Powell procedure ... \n");
  
  double tau1 = eps, tau2 = eps, eps1 = eps, eps2 = eps, const1 = 1.0e-17;
  
  double ramda = 0.05, fv, s1, s2;
  double h[8][8], s[8] = {0}, dx[8] = {0}, g0[8] = {0},
    g[8] = {0}, dg[8], wrk[8];
  
  // Initial estimate of inverse of hessian matrix
  for (int i = 0; i < 8; i++)
    for (int j = 0; j < 8; j++)
      h[i][j] = ihess(i, j);
  
  mloglikGrMP(tht, &fv, g, nthreads);
  
  if (verbose)
  {
    Rprintf("Function Value = %8.4f\n", fv);
    for (int i = 0; i < 8; ++i)
      Rprintf("Gradient[%d] = %8.2f\ttheta[%d] = %2.6f\n",
              i + 1, g[i], i + 1, tht[i]);
  }
  
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
          double sum = 0;
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
          fvout = -fv;
          aic = 2 * (fv + 8);
          if (verbose)
            Rprintf ("loglikelihood = %8.5f\tAIC = %8.5f\n", -fv, 2 * (fv + 8));
          for( int i = 0; i < 8; i++ )
          {
            dfvout[i] = g[i];
            estimate[i] = tht[i];
            for (int j = 0; j < 8; j++)
              ihess(i, j) = h[i][j];
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
          for (int i = 0; i < 8; i++)
            for (int j = i; j < 8; j++)
            {
              h[i][j] -= (dx[i] * wrk[j] + wrk[i] * dx[j] -
                dx[i] * dx[j] * (1 + s1 / s2)) / s2;
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
      double ss = 0;
      for (int i = 0; i < 8; i++)
      {
        double sum = 0;
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
      
      if ((fabs(s1) / sqrt(s2) <= tau1) && (sqrt(s2) <= tau2))
      {
        fvout = -fv;
        aic = 2 * (fv + 8);
        if (verbose)
          Rprintf ("loglikelihood = %8.5f\tAIC = %8.5f\n",
                   -fv, 2 * (fv + 8));
        for( int i = 0; i < 8; i++ )
        {
          dfvout[i] = g[i];
          estimate[i] = tht[i];
          for (int j = 0; j < 8; j++)
            ihess(i, j) = h[i][j];
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
        for (int i = 0; i < 8; i++)
        {
          for (int j = 0; j < 8; j++)
            h[i][j] = 0.0;
          h[i][i] = 1.0;
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
      for (int i = 0; i < 8; i++)
      {
        dx[i] = s[i] * ramda;
        s1 += dx[i] * dx[i];
        g0[i] = g[i];
        tht[i] += dx[i];
      }
      
      double fv0 = fv;
      mloglikGrMP(tht, &fv, g, nthreads);
      
      if (verbose)
      {
        Rprintf("Function Value = %8.4f\n", fv);
        for (int i = 0; i < 8; ++i)
          Rprintf("Gradient[%d] = %8.2f\ttheta[%d] = %2.6f\n",
                  i + 1, g[i], i + 1, tht[i]);
      }
      
      s2 = 0;
      for (int i = 0; i < 8; i++)
        s2 += g[i] * g[i];
      if (sqrt(s2) > tau2)
        continue;
      if (fv0/fv - 1 < eps1 && sqrt(s1) < eps2)
      {
        fvout = -fv;
        aic = 2 * (fv + 8);
        if (verbose)
          Rprintf ("loglikelihood = %8.5f\tAIC = %8.5f\n", -fv, 2 * (fv + 8));
        for( int i = 0; i < 8; i++ )
        {
          dfvout[i] = g[i];
          estimate[i] = tht[i];
          for (int j = 0; j < 8; j++)
            ihess(i, j) = h[i][j];
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
            int nthreads)
{
  etas data;
  data.set(revents, rpoly, tperiod, rinteg0, ndiv);
  
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

  double gparam[] = {c, p};
  double fparam[] = {D, gamma, q};

  NumericVector out(tv.length());
  double s = 0; //= mu * bk[j];
  
  for (int j = 0; j < tv.length(); j++)
  {
    int i = 0;
    while (t[i] < tv[j])
    {
      s += A * exp(alpha * m[i]) *
        gfun(tv[j] - t[i], gparam) *
        ffun(dist2(xv[j], yv[j], x[i], y[i]), m[i], fparam);
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
                int ndiv)
{
  NumericVector t = revents( _, 0), x = revents( _, 1), y = revents( _, 2),
    m = revents( _, 3), flag = revents( _, 4), bk = revents( _, 5),
    pb = revents( _, 6), lam = revents( _, 7);
  NumericVector px = rpoly( _, 0), py = rpoly( _, 1);
  
  int N = t.length(), np = px.length();
  
  // extract time period information
  const double tstart2 = tperiod[0], tlength = tperiod[1];
  
  const double mu = param[0];
  const double A = param[1];
  const double c = param[2];
  const double alpha = param[3];
  const double p = param[4];
  const double D = param[5];
  const double q= param[6];
  const double gamma = param[7];
  
  double gparam[] = {c, p};
  double fparam[] = {D, gamma, q};

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
    
    double s_thread = mu * bk[i];
    for (int j = 0; j < i; ++j)
    {
      s_thread += A * exp(alpha * m[j]) *
        gfun(t[i] - t[j], gparam) *
        ffun(dist2(x[i], y[i], x[j], y[j]), m[j], fparam);
    }
    
    lam[i] = s_thread;
    
    revents(i, 5) = bk[i];
    revents(i, 6) = mu * bk[i] / lam[i]; // probability of event i being a background event
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
              NumericVector gy)
{
  NumericVector t = revents( _, 0), x = revents( _, 1), y = revents( _, 2),
    m = revents( _, 3), pb = revents( _, 6);
  
  // extract time period information
  const double tstart2 = tperiod[0], tlength = tperiod[1];
  
  const double mu = param[0];
  const double A = param[1];
  const double c = param[2];
  const double alpha = param[3];
  const double p = param[4];
  const double D = param[5];
  const double q= param[6];
  const double gamma = param[7];
  
  double gparam[] = {c, p};
  double fparam[] = {D, gamma, q};


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
      lamb(i, j) = mu * bkgd(i, j);
      
      for (int l = 0; l < N; l++)
      {
        lamb(i, j) += A * exp(alpha * m[l]) *
          gfun(tlength - t[l], gparam) *
          ffun(dist2(x[l], y[l], gx[i], gy[j]), m[l], fparam);
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
                           int ndiv)
{
  NumericVector t = revents( _, 0), x = revents( _, 1), y = revents( _, 2),
    m = revents( _, 3);
  NumericVector px = rpoly( _, 0), py = rpoly( _, 1);
  const double tstart2 = tperiod[0], tlength = tperiod[1];
  
  const double mu = theta[0], A = theta[1], c = theta[2], alpha = theta[3];
  const double p = theta[4], D = theta[5], q = theta[6], gamma = theta[7];

  double gparam[] = {c, p};
  double fparam[] = {D, gamma, q};

  const int N = revents.nrow();
  NumericVector sinteg(N), out(N);
  
  for (int i=0; i < N; i++)
  {
    sinteg[i] =  A * exp(alpha * m[i]) *
      polyintegXX(frfunint, m[i], fparam, px, py, x[i], y[i], ndiv);
  }
  
  for (int j=0; j < N; ++j)
  {
    double sum = 0;
    for (int i=0; i < j; i++)
    {
      if (t[i] > tstart2)
        sum += gfunint(t[j] - t[i], gparam) * sinteg[i];
      else
        sum += (gfunint(t[j] - t[i], gparam) -
          gfunint(tstart2 - t[i], gparam)) * sinteg[i];
    }
    out[j] = mu * integ0 * (t[j] - tstart2) / (tlength - tstart2) + sum;
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
                           int ndiv)
{
  NumericVector t = revents( _, 0), x = revents( _, 1), y = revents( _, 2),
    m = revents( _, 3);
  NumericVector px = rpoly( _, 0), py = rpoly( _, 1);
  const double tstart2 = tperiod[0], tlength = tperiod[1];
  
  const double mu = theta[0], A = theta[1], c = theta[2], alpha = theta[3];
  const double p = theta[4], D = theta[5], q = theta[6], gamma = theta[7];
  
  double gparam[] = {c, p};
  double fparam[] = {D, gamma, q};

  const int N = revents.nrow();
  NumericVector sinteg(N);
  const int ng = tg.length();
  NumericVector out(ng);
  
  for (int i=0; i < N; i++)
  {
    sinteg[i] =  A * exp(alpha * m[i]) *
      polyintegXX(frfunint, m[i], fparam, px, py, x[i], y[i], ndiv);
  }
  
  for (int j=0; j < ng; ++j)
  {
    double sum = 0;
    for (int i=0; i < N; i++)
    {
      if (t[i] < tg[j])
      {
        sum += gfun(tg[j] - t[i], gparam) * sinteg[i];
      }
    }
    out[j] = mu * integ0 /(tlength - tstart2) + sum;
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
                          NumericVector bwd)
{
  NumericVector t = revents( _, 0), x = revents( _, 1), y = revents( _, 2),
    m = revents( _, 3), bk = revents( _, 5), pb = revents( _, 6);
  NumericVector px = rpoly( _, 0), py = rpoly( _, 1);
  const double tstart2 = tperiod[0], tlength = tperiod[1];
  
  const double mu = theta[0], A = theta[1], c = theta[2], alpha = theta[3];
  const double p = theta[4], D = theta[5], q = theta[6], gamma = theta[7];
  
  double gparam[] = {c, p};
  double fparam[] = {D, gamma, q};

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
        gint = gfunint(tlength - t[i], gparam);
      }
      else
      {
        gint = gfunint(tlength - t[i], gparam) -
          gfunint(tstart2 - t[i], gparam);
      }
      double r2 = dist2(xg[j], yg[j], x[i], y[i]);
      sum += A * exp(alpha * m[i]) * gint * ffun(r2, m[i], fparam);
      s1 += exp(-r2/(2 * bwd[i] * bwd[i])) / (2 * M_PI * bwd[i] * bwd[i]);
      s2 += pb[i] *  s1;
    }
    
    out[j] =  sum + mu * s2/(tlength - tstart2);
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
