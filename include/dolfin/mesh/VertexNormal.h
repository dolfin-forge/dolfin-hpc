// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009.
// Modified by Aurélien Larcher, 2012-16. (partial rewrite... then full rewrite)
//
// This version fixes several important issues compared to the original:
// - the basis was originally not orthogonal, not normal and randomly oriented.
// - the weights for facet normals were incorrect.
// - the maximum critical angle was hardcoded to Pi/6 so that many points were
//   set as no-slip in the case of complex geometries.
// - the original code did not consider subdomains such that normals located at
//   subdomain boundaries (like outflows) would have incorrect orientation and
//   thus incorrect enforcement of slip boundary conditions.
// - the code did not use adjacency of ranks because of two mesh distribution
//   bugs and thus sent all vertices to all ranks.
//
// The code computes normals correctly since 2013 when I reviewed it and was
// later improved.
//
// First added:  2007-05-01
// Last changed: 2015-05-27

#ifndef __DOLFIN_VERTEX_NORMAL_H
#define __DOLFIN_VERTEX_NORMAL_H

#include <dolfin/common/constants.h>
#include <dolfin/common/Array.h>
#include <dolfin/mesh/MeshFunction.h>
#include <map>

namespace dolfin
{

class BoundaryMesh;
class Mesh;
class SubDomain;
class Vertex;

/**
 *  DOCUMENTATION:
 *
 *  @class  VertexNormal
 *
 *  @brief  Provides an orthonormal basis at each vertex located on an exterior
 *          facet of the mesh, defining an outward normal vector and two
 *          tangential vectors.
 */

class VertexNormal
{

public:

  enum Type
  {
    none, unit, facet
  };

  /// Copy constructor
  VertexNormal(VertexNormal& other);

  /// Create normal, tangents for the boundary of mesh
  VertexNormal(Mesh& mesh, Type weight);

  /// Create normal, tangents for the boundary of mesh given a subdomain
  VertexNormal(Mesh& mesh, SubDomain const& subdomain, Type weight);

  /// Destructor
  ~VertexNormal();

  /// Assignment
  VertexNormal& operator=(VertexNormal& other);

  ///
  Mesh& mesh();

  ///
  Array<MeshFunction<real> *>& basis();

  ///
  MeshFunction<uint>& vertex_type();

private:

  // Cleanup
  void clear();

  // Compute normals to the boundary nodes
  void computeNormal(Mesh& mesh);

  //
  void getFacetData(VertexNormal::Type type, Mesh& mesh, BoundaryMesh& boundary,
                    Vertex& bvertex, Array<real>& normals,
                    Array<real>& weights);

  //--- ATTRIBUTES ------------------------------------------------------------

  // Global mesh
  Mesh& mesh_;

  SubDomain const * const subdomain_;

  //
  Array<MeshFunction<real> *> basis_;

  // Define vertex type: 1 surface, 2 edge, 3 surface
  MeshFunction<uint> vertex_type_;

  // Maximum absolute angle between two neighbouring facets
  real const alpha_max_;

  Type type_;

};

//-----------------------------------------------------------------------------
inline Mesh& VertexNormal::mesh()
{
  return mesh_;
}

//-----------------------------------------------------------------------------
inline Array<MeshFunction<real> *>& VertexNormal::basis()
{
  return basis_;
}

//-----------------------------------------------------------------------------
inline MeshFunction<uint>& VertexNormal::vertex_type()
{
  return vertex_type_;
}

}
#endif
