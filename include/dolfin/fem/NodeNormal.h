// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009.
// Modified by Aurélien Larcher, 2012-13. (partial rewrite)
//
// First added:  2007-05-01
// Last changed: 2009-03-17

#ifndef __NODENORMAL_H
#define __NODENORMAL_H

#include <dolfin/fem/BoundaryNormal.h>
#include <dolfin/mesh/Mesh.h>

#include <dolfin/common/constants.h>
#include <dolfin/common/Array.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <map>

namespace dolfin
{
class NodeNormal : public BoundaryNormal
{
public:

  /// Copy constructor
  NodeNormal(NodeNormal& node_normal);

  /// Create normal, tangents for the boundary of mesh
  NodeNormal(Mesh& mesh);

  ///
  ~NodeNormal();

  ///
  void compute();

  /// Assignment
  NodeNormal& operator=(NodeNormal& node_normal);

  /// Define mesh functions for normal and tangents
  /// These are merely aliases now
  MeshFunction<real> * normal;
  MeshFunction<real> * tau;
  MeshFunction<real> * tau_1;
  MeshFunction<real> * tau_2;

  /// Define vertex type as the number of discriminated hyperplanes:
  /// 1 surface, 2 edge, >= 3 corner
  MeshFunction<uint> vertex_type;

private:

  /// Cleanup
  void Clear();

  /// Compute normals to the boundary nodes
  void ComputeNormal(Mesh& mesh);

  ///
  void ComputeTangentialVectors(Mesh& mesh, Function& Fnormal, Function& Ftau,
                                NodeNormal& node_normal);

  ///
  void ComputeTangentialVectors(Mesh& mesh, Function& Fnormal, Function& Ftau_1,
                                Function& Ftau_2, NodeNormal& node_normal);

  ///
  void CacheSharedArea(Mesh& mesh, BoundaryMesh& boundary);

  //--- ATTRIBUTES ------------------------------------------------------------

  Mesh& mesh;

  Array<MeshFunction<real> *> basis_;

  /// Entities shared between processors
  Array<real> shared_vertexnormals_;
  std::map<uint, Array<real> > shared_facetnormals_block_;
  std::map<uint, Array<real> > shared_facetweights_block_;

  /// Number of boundary mesh cells (facets for global) neighbouring a boundary
  /// vertex
  std::map<uint, uint> num_neigh_cells_;
  std::map<uint, uint> shared_offsetidx_;
  uint vertex_offset_;
  uint facetnormals_offset_;
  uint facetweights_offset_;

  /// Should be set to the size of the offset information stored for each vertex
  /// Padding = 3: (NbNeighbouringCells, FacetNormalOffset, FacetWeightOffset)
  static uint const offsetidx_padding_ = 3;

  /// Maximum absolute angle between two neighbouring facets to be discriminated
  /// as belonging to difference hyperplanes.
  real const alpha_max_;

  /// Weighing use to computing the node normal from facet normals.
  enum weight_type
  {
    none, facet, cell
  };

  weight_type weighting_;

};

}
#endif

