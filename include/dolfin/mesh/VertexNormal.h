// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009.
// Modified by Aurélien Larcher, 2012-13. (partial rewrite)
//
// First added:  2007-05-01
// Last changed: 2009-03-17

#ifndef __VERTEXNORMAL_H
#define __VERTEXNORMAL_H

#include <dolfin/common/constants.h>
#include <dolfin/common/Array.h>
#include <dolfin/mesh/MeshFunction.h>
#include <map>

namespace dolfin
{

class BoundaryMesh;
class Mesh;

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
    none, facet, cell
  };

  /// Copy constructor
  VertexNormal(VertexNormal& other);

  /// Create normal, tangents for the boundary of mesh
  VertexNormal(Mesh& mesh, Type weight);

  /// Destructor
  ~VertexNormal();

  /// Assignment
  VertexNormal& operator=(VertexNormal& other);

  ///
  Mesh& mesh();

  ///
  Array<MeshFunction<real> *> const& basis() const;

  ///
  MeshFunction<uint> const& vertex_type() const;

private:

  // Cleanup
  void Clear();

  // Compute normals to the boundary nodes
  void ComputeSimpleNormal(Mesh& mesh);

  // Compute normals to the boundary nodes
  void ComputeNormal(Mesh& mesh);

  ///
  void CacheSharedArea(Mesh& mesh, BoundaryMesh& boundary);

  /// Implemented as a 3D vector.
  void NormalizeVector(real (&v)[3]);

  ///
  void GetLocalFacetsData(uint const& gdim, Vertex& vertex,
                          MeshFunction<uint>& cell_map, uint& nb_neigh,
                          Array<real>& normals, Array<real>& weights);

  //--- ATTRIBUTES ------------------------------------------------------------

  // Global mesh
  Mesh& mesh_;

  //
  Array<MeshFunction<real> *> basis_;

  // Define vertex type: 1 surface, 2 edge, 3 surface
  MeshFunction<uint> vertex_type_;

  //
  Array<real> shared_normal;
  std::map<uint, Array<real> > shared_facetnormals_block_;
  std::map<uint, Array<real> > shared_facetweights_block_;

  // Number of boundary mesh cells (facets for global) neighbouring a boundary
  // vertex
  std::map<uint, uint> num_neigh_cells_;
  std::map<uint, uint> shared_offsetidx_;
  uint vertex_offset_;
  uint facetnormals_offset_;
  uint facetweights_offset_;

  // Should be set to the size of the offset information stored for each vertex
  // Padding = 3: (NbNeighbouringCells, FacetNormalOffset, FacetWeightOffset)
  static uint const offsetidx_padding_ = 3;

  // Maximum absolute angle between two neighbouring facets
  real const alpha_max_;

  Type weighting_;

};

//-----------------------------------------------------------------------------
inline Mesh& VertexNormal::mesh()
{
  return mesh_;
}

//-----------------------------------------------------------------------------
inline Array<MeshFunction<real> *> const& VertexNormal::basis() const
{
  return basis_;
}

//-----------------------------------------------------------------------------
inline MeshFunction<uint> const& VertexNormal::vertex_type() const
{
  return vertex_type_;
}

//-----------------------------------------------------------------------------
inline void VertexNormal::NormalizeVector(real (&v)[3])
{
  real nrm = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  dolfin_assert(nrm > 0);
  v[0] /= nrm;
  v[1] /= nrm;
  v[2] /= nrm;
}

}
#endif
