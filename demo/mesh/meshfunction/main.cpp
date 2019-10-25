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
  // Read mesh from file
  Mesh mesh("mesh2D.xml.gz");

  // Read mesh function from file
  MeshValues<real, Cell> mv(mesh);
  File f("meshfunction.xml");
  f >> mv;
  // MeshFunction<real> f(mesh, "meshfunction.xml");

  // Write mesh function to file
  File out("meshfunction_out.pvd");
  out << mv;
}
