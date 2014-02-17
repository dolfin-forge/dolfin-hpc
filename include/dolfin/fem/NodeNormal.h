// Copyright (C) 2014 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-30
// Last changed: 2014-01-30

#ifndef __NODENORMAL_H
#define __NODENORMAL_H

#include <dolfin/fem/BoundaryNormal.h>

#include <dolfin/common/constants.h>
#include <dolfin/common/Array.h>
#include <dolfin/function/Function.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/VertexNormal.h>

namespace dolfin
{

class FiniteElementSpace;
class Mesh;

/**
 *  DOCUMENTATION:
 *
 *  @class  NodeNormal
 *
 *  @brief  Provides an orthonormal basis at each geometrical node located on an
 *          exterior facet of the mesh for a given finite element space,
 *          defining an outward normal vector and two tangential vectors.
 */

class NodeNormal : public BoundaryNormal
{
public:

  /// Create normal, tangents to the boundary of mesh at vertices
  NodeNormal(Mesh& function, VertexNormal::Type weight);

  /// Destructor
  ~NodeNormal();

  ///
  void init(FiniteElementSpace& space);

private:

  /// Cleanup
  void Clear();

  /// Compute normals to the boundary nodes
  void ComputeBasisP1();

  //--- ATTRIBUTES ------------------------------------------------------------

  uint const tdim_;
  mutable FiniteElementSpace const * space_;
  mutable bool local_space_;
  VertexNormal normals_;
  Array<MeshFunction<real> *> const& meshbasis_;
  Array<GenericVector *> V_;

};

}
#endif

