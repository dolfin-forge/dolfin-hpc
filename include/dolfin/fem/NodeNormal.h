// Copyright (C) 2014 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-30
// Last changed: 2014-01-30

#ifndef __NODENORMAL_H
#define __NODENORMAL_H

#include <dolfin/common/constants.h>
#include <dolfin/common/Array.h>
#include <dolfin/function/Function.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/VertexNormal.h>

namespace dolfin
{

class BoundaryMesh;
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

class NodeNormal
{
public:

  /// Create normal, tangents to the boundary of mesh at vertices
  NodeNormal(Function& function, VertexNormal::Type weight);

  /// Destructor
  ~NodeNormal();

  /// Return mesh
  Mesh& mesh()
  {
    return mesh_;
  }

  /// Return the orthonormal basis (n, tau) in 2d or (n, tau1, tau2) in 3d
  Array<Function> const& basis() const;

private:

  /// Cleanup
  void Clear();

  /// Compute normals to the boundary nodes
  void ComputeBasisP1();

  //--- ATTRIBUTES ------------------------------------------------------------

  Mesh& mesh_;
  uint const tdim_;
  mutable FiniteElementSpace const * space_;
  mutable bool local_space_;
  VertexNormal normals_;
  Array<MeshFunction<real> *> const& meshbasis_;
  Array<Function> basis_;
  Array<GenericVector *> V_;
  Function nodetype_;

};

//-----------------------------------------------------------------------------
inline Array<Function> const& NodeNormal::basis() const
{
  return basis_;
}

}
#endif

