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

#ifdef ENABLE_UFL
#include "ufc2/Lift.h"
#include "ufc2/Drag.h"
#else
#include "ufc1/Lift.h"
#include "ufc1/Drag.h"
#endif

using namespace dolfin;

int main()
{
  // Read velocity field from file and get the mesh
#ifdef ENABLE_UFL 
  Function p("ufc2/pressure.xml.gz");
  Mesh& mesh(p.mesh());
#else
  Function p("ufc1/pressure.xml.gz");
  Mesh& mesh(p.mesh());
#endif

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
  
  // Facet normal
//  FacetNormal n(mesh);

  // Functionals for lift and drag
  std::map<std::string const, Function *> coef_map;
  coef_map["p"] = &p;
  LiftFunctional L(mesh, coef_map);
  DragFunctional D(mesh, coef_map);

  // Assemble functionals over sub domain
  Fish fish;
  Assembler assembler(mesh);
  real lift = assembler.assemble(L, fish);
  real drag = assembler.assemble(D, fish);

  message("Lift: %f", lift);
  message("Drag: %f", drag);
}
