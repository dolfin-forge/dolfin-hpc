// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-03-11
// Last changed: 2007-03-17
//
// Testing evaluation at arbitrary points

#include <dolfin.h>
#include "Projection.h"

using namespace dolfin;

class F : public ScalarExpression
{
public:
  
  F() : ScalarExpression() {}

  void eval(real* values, const real* x) const
  {
    values[0] = sin(3.0*x[0])*sin(3.0*x[1])*sin(3.0*x[2]);
  }
};

int main()
{
  // Create mesh and a point in the mesh
  UnitCube mesh(8, 8, 8);
  real x[3] = {0.3, 0.3, 0.3};
  real f_values[1] = {0.0};
  real g_values[1] = {0.0};

  // A user-defined function
  F force;
  Function f(mesh, force);

  // Project to a discrete function
  ProjectionBilinearForm a(mesh);
  ProjectionLinearForm L(f);
  LinearPDE pde(a, L, mesh);
  Function g(mesh);
  pde.solve(g);

  // Evaluate user-defined function f
  f.eval(f_values, x);
  message("f(x) = %g", f_values[0]);

  // Evaluate discrete function g (projection of f)
  g.eval(g_values, x);
  message("g(x) = %g", g_values[0]);
}
