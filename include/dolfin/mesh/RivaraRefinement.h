// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009-2010.
//

#ifndef __DOLFIN_RIVARA_REFINEMENT_H
#define __DOLFIN_RIVARA_REFINEMENT_H

#include <dolfin/common/types.h>
#include <dolfin/mesh/MeshFunction.h>

#include <dolfin/config/dolfin_config.h>

namespace dolfin
{
  class Mesh;

  class RivaraRefinement
  {
  public:
    
    /// Refine simplicial mesh locally by recursive edge bisection 
    static void refine(Mesh& mesh, 
                       MeshFunction<bool>& cell_marker,
                       real tf = 0.0, 
                       real tb = 0.0, 
                       real ts = 0.0,
                       bool balance = true);
  };
}
#endif
