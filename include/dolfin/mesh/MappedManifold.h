// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// This class is the first implementation allowing identification of mesh
// boundaries through periodic mappings on general meshes.
//
// First added:  2014-10-05
// Last changed: 2014-10-05

#ifndef __DOLFIN_MAPPED_MANIFOLD_H
#define __DOLFIN_MAPPED_MANIFOLD_H

#include <dolfin/common/types.h>

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshDependent.h>

namespace dolfin
{

class Cell;
class PeriodicSubDomain;
class Vertex;

/**
 *  @class  MappedManifold
 */

class MappedManifold : public Mesh, public MeshDependent
{

public:

  /// Create boundary mesh from given mesh
  MappedManifold(Mesh& mesh, PeriodicSubDomain const& subdomain);

  /// Destructor
  ~MappedManifold();

  /// Return facet index in the mesh associated with the boundary cell
  uint facet_index(Cell const& boundary_cell);

  /// Return vertex index in the mesh associated with the boundary vertex
  uint vertex_index(Vertex const& boundary_vertex);

  /// Return periodic domain that generated the instance
  PeriodicSubDomain const& subdomain() const
  {
    return subdomainG_;
  }

  /// Return facets contained in G
  _set<uint> const& Gfacets() const
  {
    return facetsG_;
  }

  /// Return facets contained in H
  _set<uint> const& Hfacets() const
  {
    return facetsH_;
  }

  /// Return facets contained in G inter H
  _set<uint> const& Ifacets() const
  {
    return facetsI_;
  }

  /// Return H facets for which image overlaps with one local G facet
  _set<uint> const& Lfacets() const
  {
    return facetsL_;
  }

private:

  void init();

  PeriodicSubDomain const& subdomainG_;

  _set<uint> facetsG_;
  _set<uint> facetsH_;
  _set<uint> facetsI_; // Intersection of G and H
  _set<uint> facetsL_; // H facet with local G facet

  Array<uint> vertex_map_;
  Array<uint> cell_map_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_MAPPED_MANIFOLD_H */
