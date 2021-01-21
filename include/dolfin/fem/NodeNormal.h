// Copyright (C) 2014 Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.

// The original class named NodeNormal from UNICORN is actually VertexNormal
// since the former was written for linear Lagrange in 3D only to compute a
// non-orthogonal, non-anticlockwise, sometimes non-normal basis with no
// support for subdomains.
// This class together with the rewrite of the finite element framework
// provides support for the computation of a normal field and tangential vectors
// interpolated to a Lagrange finite element space of arbitrary order.

#ifndef __DOLFIN_NODENORMAL_H
#define __DOLFIN_NODENORMAL_H

#include <dolfin/common/Array.h>
#include <dolfin/common/constants.h>
#include <dolfin/fem/BoundaryNormal.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/parameter/parameters.h>

namespace dolfin
{

class SubDomain;

class NodeNormal : public BoundaryNormal
{

public:

  enum Type
  {
    none, unit, facet
  };

  /// Create normal, tangents for the boundary of mesh
  NodeNormal( Mesh& mesh, Type w = unit,
              real alpha = dolfin_get<real>("NodeNormal alpha") );

  /// Create normal, tangents for the boundary of mesh for given subdomain
  NodeNormal( Mesh& mesh, SubDomain const& subdomain, Type w = none,
              real alpha = dolfin_get<real>("NodeNormal alpha") );

  /// Destructor
  ~NodeNormal() override = default;

  /// Compute the orthogonal basis
  void compute() override;

  /// Returns the node type
  auto node_type(uint node_id) const -> uint;

private:

  /// Cleanup
  void clear();

  /// Compute boundary normal basis
  void compute(Mesh& mesh, Array<Function>& basis);

  //--- ATTRIBUTES ------------------------------------------------------------

  Mesh& mesh_;

  SubDomain const * const subdomain_;

  /// Maximum absolute angle between two neighbouring facets to be discriminated
  /// as belonging to different hyperplanes.
  real const alpha_max_;

  /// Type of weight used for computing the node normal from facet normals.
  Type const type_;

  ///
  _map<uint, uint> node_type_;
};

}
#endif

