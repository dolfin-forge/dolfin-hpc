// Copyright (C) 2006 Ola Skavhaug.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2007.
//
// First added:  2006-11-29
// Last changed: 2008-03-31

#include <dolfin.h>

using namespace dolfin;

int main()
{
  dolfin_init();

  // Read mesh from file
  Mesh mesh("mesh2D.bin");

  // Read mesh function from file
  MeshValues<real, Cell> mv(mesh);
  File("meshfunction.bin") >> mv;

  // Write mesh function to file
  File("meshfunction_out.pvd") << mv;

  dolfin_finalize();

  return 0;
}
