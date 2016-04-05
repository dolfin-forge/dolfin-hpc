// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009.
// Modified by Aurélien Larcher, 2012-15. (partial rewrite... then rewrite)
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
 *  @class  NodeNormal
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

  ///
  static uint computeBasis(uint gdim, Point B[], Array<real> N, Array<real> W,
                           real cosalpha_max, bool weighted);

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

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_VERTEX_NORMAL_H */
