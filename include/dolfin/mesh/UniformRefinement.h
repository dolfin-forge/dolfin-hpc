// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-06-07
// Last changed: 2006-06-16

#ifndef __DOLFIN_MESH_UNIFORM_REFINEMENT_H
#define __DOLFIN_MESH_UNIFORM_REFINEMENT_H

#include <dolfin/common/types.h>

namespace dolfin
{

class Mesh;

/// This class implements uniform mesh refinement for different mesh types.

struct UniformRefinement
{

  /// Refine mesh uniformly according to default pattern
  void operator()(Mesh& mesh);

};

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_UNIFORM_REFINEMENT_H */
