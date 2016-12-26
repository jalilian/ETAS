#include <Rinternals.h>

// *******************************************************************************

double fr(double r, double w[]);

double dgamma_fr(double r, double w[]);

double dD_fr(double r, double w[]);

double dq_fr(double r, double w[]);

// *******************************************************************************

double clambdaj(double *theta,
               int j,
               double *t,
               double *x,
               double *y,
               double *m,
               double *bk);

void clambdajGr(double *theta,
               int j,
               double *t,
               double *x,
               double *y,
               double *m,
               double *bk,
	       double *fv,
	       double *dfv);

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
              double *tlength);

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
	      double *dfv);

// *******************************************************************************

