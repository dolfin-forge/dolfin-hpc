// Copyright (C) 2004-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells 2006, 2007.
// Modified by Dag Lindbo, 2008.
// Modified by Aurélien Larcher, 2014.
//
// First added:  2004
// Last changed: 2007-12-09

#ifndef __DOLFIN_LINEAR_PDE_H
#define __DOLFIN_LINEAR_PDE_H

#include <dolfin/parameter/Parametrized.h>

#include <dolfin/common/Array.h>

namespace dolfin
{

class BilinearForm;
class LinearForm;
class Form;
class Mesh;
class BoundaryCondition;
class Function;

/**
 *  @class  LinearPDE
 *
 *  @brief  A LinearPDE represents a (system of) linear partial differential
 *          equation(s) in variational form: Find u in V such that
 *
 *               a(v, u) = L(v) for all v in V',
 *
 *          where a is a bilinear form and L is a linear form.
 */

class LinearPDE : public Parametrized
{
public:

  /// Define a linear PDE with natural boundary conditions
  LinearPDE(BilinearForm& a, LinearForm& L);

  /// Define a linear PDE with a single Dirichlet boundary condition
  LinearPDE(BilinearForm& a, LinearForm& L, BoundaryCondition& bc);

  /// Define a linear PDE with a set of Dirichlet boundary conditions
  LinearPDE(BilinearForm& a, LinearForm& L, Array<BoundaryCondition*>& bcs);

  /// Destructor
  ~LinearPDE();

  /// Solve PDE system
  void solve(Function& u);

private:

  // The bilinear form
  BilinearForm& a;

  // The linear form
  LinearForm& L;

  // The boundary conditions
  Array<BoundaryCondition*> bcs;

  // Not assembled ?
  bool not_assembled;
};

}

#endif
