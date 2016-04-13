// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-06-07
// Last changed: 2006-06-16

#ifndef __DOLFIN_UNIFORM_MESH_REFINEMENT_H
#define __DOLFIN_UNIFORM_MESH_REFINEMENT_H

#include <dolfin/common/types.h>

namespace dolfin
{

class Mesh;

/// This class implements uniform mesh refinement for different mesh types.

class UniformMeshRefinement
{

public:

  /// Refine mesh uniformly according to mesh type
  static void refine(Mesh& mesh);

private:

  /// Refine simplicial mesh uniformly
  static void refineSimplex(Mesh& mesh);

  /// Refine hypercubic mesh uniformly
  static void refineHypercube(Mesh& mesh);

};

} /* namespace dolfin */

#endif /* __DOLFIN_UNIFORM_MESH_REFINEMENT_H */
