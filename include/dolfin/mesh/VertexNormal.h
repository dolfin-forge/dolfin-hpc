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

  // Copy constructor
  VertexNormal(VertexNormal& other);

  // Create normal, tangents for the boundary of mesh
  VertexNormal(Mesh& mesh);

  /// Destructor
  ~VertexNormal();

  // Assignment
  VertexNormal& operator=(VertexNormal& other);

  // Define mesh functions for normal and tangents
  // These are merely aliases now
  // One-liner is bad coding style.
  MeshFunction<real> * normal;
  MeshFunction<real> * tau;
  MeshFunction<real> * tau_1;
  MeshFunction<real> * tau_2;

  // Define vertex type: 1 surface, 2 edge, 3 surface
  MeshFunction<uint> vertex_type;

  Array<MeshFunction<real> *> const& basis() const;

private:

  // Cleanup
  void Clear();

  // Compute normals to the boundary nodes
  void ComputeNormal(Mesh& mesh);

  ///
  void CacheSharedArea(Mesh& mesh, BoundaryMesh& boundary);

  //--- ATTRIBUTES ------------------------------------------------------------

  Mesh& mesh;

  Array<MeshFunction<real> *> basis_;

  Array<real> shared_vertexnormals_;
  std::map<uint, Array<real> > shared_facetnormals_block_;
  std::map<uint, Array<real> > shared_facetweights_block_;

  // Number of boundary mesh cells (facets for global) neighbouring a boundary vertex
  std::map<uint, uint> num_neigh_cells_;
  std::map<uint, uint> shared_offsetidx_;
  uint vertex_offset_;
  uint facetnormals_offset_;
  uint facetweights_offset_;

  // Should be set to the size of the offset information stored for each vertex
  // Here padding = 3: (NbNeighbouringCells, FacetNormalOffset, FacetWeightOffset)
  static uint const offsetidx_padding_ = 3;

  // Maximum absolute angle between two neighbouring facets
  real const alpha_max_;

  enum weight_type
  {
    none, facet, cell
  };

  weight_type weighting_;

};

//-----------------------------------------------------------------------------
inline Array<MeshFunction<real> *> const& VertexNormal::basis() const
{
  return basis_;
}

}
#endif
