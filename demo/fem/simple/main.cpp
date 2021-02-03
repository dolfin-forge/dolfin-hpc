// Copyright (C) 2003-2005 Johan Hoffman and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2007.
//
// Thanks to David Heintz for the reference matrices.
//
// This demo program dmonstrates how to create simple finite
// element matrices like the stiffness matrix and mass matrix.
// For general forms and matrices, forms must be defined and
// compiled with FFC.

#include <dolfin.h>

using namespace dolfin;

int main()
{
  dolfin_init();

  // Load reference mesh (just a simple tetrahedron)
  Mesh mesh("tetrahedron.bin");

  // Create stiffness and mass matrices
  // StiffnessMatrix A(mesh);
  // MassMatrix M(mesh);

  // Create reference matrices
  real A0_array[4][4] = { {  1.0/2.0, -1.0/6.0, -1.0/6.0, -1.0/6.0 },
                          { -1.0/6.0,  1.0/6.0,  0.0,      0.0     },
                          { -1.0/6.0,  0.0,      1.0/6.0,  0.0     },
                          { -1.0/6.0,  0.0,      0.0,      1.0/6.0 } };

  real M0_array[4][4] = { { 1.0/60.0,  1.0/120.0, 1.0/120.0, 1.0/120.0 },
                          { 1.0/120.0, 1.0/60.0,  1.0/120.0, 1.0/120.0 },
                          { 1.0/120.0, 1.0/120.0, 1.0/60.0,  1.0/120.0 },
                          { 1.0/120.0, 1.0/120.0, 1.0/120.0, 1.0/60.0  } };

  size_t position[4] = {0, 1, 2, 3};

  Matrix A0(4,4);
  Matrix M0(4,4);
  A0.set(*A0_array, 4, position, 4, position);
  M0.set(*M0_array, 4, position, 4, position);

  A0.apply();
  M0.apply();

  std::cout << "Reference stiffness matrix:\n";
  A0.disp();

  std::cout << "Reference mass matrix:\n";
  M0.disp();

  dolfin_finalize();

  return 0;
}
