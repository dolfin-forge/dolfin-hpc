// Copyright (C) 2003-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2003-06-03
// Last changed: 2007-07-18

#include <dolfin/config/dolfin_config.h>


#include <cmath>
#include <dolfin/common/constants.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/math/Legendre.h>
#include <dolfin/quadrature/GaussianQuadrature.h>

#if HAVE_F77_LAPACK
extern "C" {
  void dgesv_(unsigned int *n, int *nrhs, double *A, unsigned int *lda, int *ipiv, 
	      double *b, unsigned int *ldb, int* info);
}
#endif
#define RM(row,col,nrow) ((row) + ((nrow)*(col)))


using namespace dolfin;

//-----------------------------------------------------------------------------
GaussianQuadrature::GaussianQuadrature(unsigned int n) : Quadrature(n)
{
  // Length of interval [-1,1]
  m = 2.0;
}
//-----------------------------------------------------------------------------
void GaussianQuadrature::init()
{
  computePoints();
  computeWeights();
}
//-----------------------------------------------------------------------------
void GaussianQuadrature::computeWeights()
{
  // Compute the quadrature weights by solving a linear system of equations
  // for exact integration of polynomials. We compute the integrals over
  // [-1,1] of the Legendre polynomials of degree <= n - 1; These integrals
  // are all zero, except for the integral of P0 which is 2.
  //
  // This requires that the n-point quadrature rule is exact at least for
  // polynomials of degree n-1.

  // Special case n = 0
  if ( n == 0 )
  {
    weights[0] = 2.0;
    return;
  }
#if HAVE_F77_LAPACK
  int nrhs = 1;
  int info;

  real *A = new real[n * n];
  real *x = new real[n];
  real *b = new real[n];
  int *ipiv = new int[n];
  
  // Compute the matrix coefficients
  for (unsigned int i = 0; i < n; i++)
  {
    Legendre p(i);
    for (unsigned int j = 0; j < n; j++)
      A[RM(i, j, n)] = p(points[j]);
    b[i] = 0.0;
  }
  b[0] = 2.0;
  
  // Solve the system of equations
  dgesv_(&n, &nrhs, &A[0], &n, &ipiv[0], &b[0], &n, &info);
  
  // Save the weights
  for (unsigned int i = 0; i < n; i++)
    weights[i] = x[i];

  delete[] ipiv;
  delete[] b;
  delete[] x;
  delete[] A;
#else
  error("GaussianQuadrature needs F77 LAPACK");
#endif
}
//-----------------------------------------------------------------------------
bool GaussianQuadrature::check(unsigned int q) const
{
  // Checks that the points and weights are correct. We compute the
  // value of the integral of the Legendre polynomial of degree q.
  // This value should be zero for q > 0 and 2 for q = 0
  
  Legendre p(q);
  
  real sum = 0.0;
  for (unsigned int i = 0; i < n; i++)
    sum += weights[i] * p(points[i]);
  
  //message("Checking quadrature weights: %.2e.", fabs(sum));
  
  if ( q == 0 )
  {
    if ( fabs(sum - 2.0) < 100.0*DOLFIN_EPS )
      return true;
  }
  else
  {
    if ( fabs(sum) < 100.0*DOLFIN_EPS )
      return true;
  }

  message("Quadrature check failed: r = %.2e.", fabs(sum));

  return false;
}
//-----------------------------------------------------------------------------


