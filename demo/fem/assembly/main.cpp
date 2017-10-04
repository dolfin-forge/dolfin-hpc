// Copyright (C) 2008 Anders Logg and Magnus Vikstrøm
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-02-04
// Last changed: 2017-10-04
//

#include <dolfin/common/Test.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/mesh/UnitCube.h>

#include "ReactionDiffusion.h"

using namespace dolfin;

int main()
{
  Test T;

  // Create mesh and form
  UnitCube mesh(20, 20, 20);
  ReactionDiffusion::BilinearForm a(mesh); 

  // Assemble matrix
  Matrix A;
  a.assemble(A, true);
  message("|A| = %e", A.norm());
 
  return 0;
}
