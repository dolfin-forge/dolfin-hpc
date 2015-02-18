// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-09-19
// Last changed: 2007-04-30
//
// This demo program computes the value of the functional
//
//     M(v) = int v^2 + (grad v)^2 dx
//
// on the unit square for v = sin(x) + cos(y). The exact
// value of the functional is M(v) = 2 + 2*sin(1)*(1-cos(1))
//
// The functional M corresponds to the energy norm for a
// simple reaction-diffusion equation.

#include <dolfin.h>

#include "EnergyNorm.h"

using namespace dolfin;

int main()
{
  // The function v
  class MyFunction : public ScalarExpression
  {
  public:

    MyFunction() : ScalarExpression() {}
    
    void eval(real* values, const real* x) const
    {
      values[0] = sin(x[0]) + cos(x[1]);
    }
  };

  // Compute approximate value
  UnitSquare mesh(16, 16);
  MyFunction v;
  Function my_function(mesh, v);
  EnergyNormFunctional M(my_function);
  Assembler assembler(mesh);
  real value = assembler.assemble(M);

  // Compute exact value
  real exact_value = 2.0 + 2.0*sin(1.0)*(1.0 - cos(1.0));

  message("The energy norm of v is %.15g (should be %.15g).", value, exact_value);
  
  return 0;
}
