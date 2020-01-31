// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-05-14
// Last changed: 2007-08-20
//
// This demo demonstrates how to compute functionals (or forms
// in general) over subsets of the mesh. The two functionals
// lift and drag are computed for the pressure field around
// a dolphin. Here, we use the pressure field obtained from
// solving the Stokes equations (see demo program in the
// sub directory src/demo/pde/stokes/taylor-hood).

#include <dolfin.h>

#include "Lift.h"
#include "Drag.h"

using namespace dolfin;

// Define sub domain for the dolphin
class Fish : public SubDomain
{
  bool inside(const real* x, bool on_boundary) const
  {
    return (x[0] > DOLFIN_EPS && x[0] < (1.0 - DOLFIN_EPS) &&
            x[1] > DOLFIN_EPS && x[1] < (1.0 - DOLFIN_EPS) &&
            on_boundary);
  }
};

int main()
{
  dolfin_init();

  // Read velocity field from file and get the mesh
  Function p;
  File("pressure.xml.gz") >> p;
  Mesh& mesh(p.mesh());

  // Functionals for lift and drag
  Lift::Functional L(mesh, p);
  Drag::Functional D(mesh, p);

  // Assemble functionals over sub domain
  Fish fish;
  Scalar l, d;
  Assembler::assemble( l, L, fish, true );
  real lift = l;

  Assembler::assemble( d, D, fish, true );
  real drag = d;


  message("Lift: %f", lift);
  message("Drag: %f", drag);

  dolfin_finalize();
}
