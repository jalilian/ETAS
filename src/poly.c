#include<Rmath.h>
#include "dist.h"

// *******************************************************************************

double frint(double (*func)(double, double []), 
	     double funcpara[], 
             double x1, 
	     double y1, 
	     double x2, 
	     double y2, 
	     double cx, 
	     double cy)
{
  int id = 1;
  double det, r0, r1, r2, r12, theta, f1, f2, f3, x0, y0;
  
  det = (x1 * y2 + y1 * cx + x2 * cy) - (x2 * y1 + y2 * cx + x1 * cy);
  
  if (det < 0)
    id = -1;
  if (fabs(det) < 1.0e-10)
    return 0;
  
  r1 = dist(x1, y1, cx, cy);
  r2 = dist(x2, y2, cx, cy);
  r12 = dist(x1, y1, x2, y2);
  theta = (r1 * r1 + r2 * r2 - r12 * r12)/(2 * r1 * r2);
  if (fabs(theta) > 1)
    theta = 1 - 1.0e-10;

  theta = acos(theta);

  if (r1 + r2 > 1.0e-20)
    {
      x0 = x1 + r1/(r1 + r2) * (x2 - x1);
      y0 = y1 + r1/(r1 + r2) * (y2 - y1);
    }
  else
    return 0;
  
  r0 = dist(x0, y0, cx, cy);
  
  f1 = func(r1, funcpara);
  f2 = func(r0, funcpara);
  f3 = func(r2, funcpara);
  
  return id * (f1/6 + (f2 * 2)/3 + f3/6) * theta;
}

// *******************************************************************************
// approximating the integral of a function on a polygon region

double polyinteg(double (*func)(double, double []), 
		 double funcpara[], 
                 int *np, 
		 double *px, 
		 double *py, 
		 double cx, 
		 double cy)
{
  int ndiv = 1000; 
  double sum = 0, dxx, dyy, x1, x2, y1, y2;
  
  for (int j = 0; j < (*np - 1); j++)
    {
      dxx = (px[j + 1] - px[j]) / ndiv;
      dyy = (py[j + 1] - py[j]) / ndiv;
      for( int i = 0; i < ndiv; i++)
	{
	  x1 = px[j] + dxx * i;
	  y1 = py[j] + dyy * i;
	  x2 = px[j] + dxx * (i + 1);
	  y2 = py[j] + dyy * (i + 1);
	  
	  sum += frint(func, funcpara, x1, y1, x2, y2, cx, cy);
	}
    }

  return sum;
}
