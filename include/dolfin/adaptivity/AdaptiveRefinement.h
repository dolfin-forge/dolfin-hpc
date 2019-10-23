// Copyright (C) 2010 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_ADAPTIVE_REFINEMENT_H
#define __DOLFIN_ADAPTIVE_REFINEMENT_H

#include <dolfin/common/types.h>
#include <dolfin/mesh/MeshValues.h>

#include <vector>

namespace dolfin
{

template<class T> class Array;
class Function;
class Mesh;
class Function;
class Vector;

class AdaptiveRefinement
{
public:

  /// Refine mesh using "simple" of "rivara" strategy
  static void refine(Mesh& mesh, MeshValues<bool, Cell>& cell_marker);

  ///
  static void refine_and_project(Mesh& mesh, Array<Function *> const& functions,
                                 MeshValues<bool, Cell>& cell_marker);

private:

  ///
  static void redistribute_func(Mesh& mesh, Function const& f, real **vp,
                                uint **rp, uint& m,
                                MeshValues<uint, Cell>& distribution);

  /// Project function on new mesh i.e. interpolation on non-matching meshes.
  static void project(Mesh& new_mesh, Array<Function *>& f_post,
                      Function& projected);

};
}
#endif
