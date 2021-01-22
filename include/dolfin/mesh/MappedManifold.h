// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

// This class is the first implementation allowing identification of mesh
// boundaries through periodic mappings on general meshes.

#ifndef __DOLFIN_MAPPED_MANIFOLD_H
#define __DOLFIN_MAPPED_MANIFOLD_H

#include <dolfin/common/types.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshDependent.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/entities/Vertex.h>

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
  ~MappedManifold() override = default;

  /// Return facet index in the mesh associated with the boundary cell
  auto facet_index(Cell const& boundary_cell) -> size_t;

  /// Return vertex index in the mesh associated with the boundary vertex
  auto vertex_index(Vertex const& boundary_vertex) -> size_t;

  /// Return periodic domain that generated the instance
  auto subdomain() const -> PeriodicSubDomain const&
  {
    return subdomainG_;
  }

  /// Return facets contained in G
  auto Gfacets() const -> _set<size_t> const&
  {
    return facetsG_;
  }

  /// Return facets contained in H
  auto Hfacets() const -> _set<size_t> const&
  {
    return facetsH_;
  }

  /// Return facets contained in G inter H
  auto Ifacets() const -> _set<size_t> const&
  {
    return facetsI_;
  }

  /// Return H facets for which image overlaps with one local G facet
  auto Lfacets() const -> _set<size_t> const&
  {
    return facetsL_;
  }

private:

  void init();

  PeriodicSubDomain const& subdomainG_;

  _set<size_t> facetsG_;
  _set<size_t> facetsH_;
  _set<size_t> facetsI_; // Intersection of G and H
  _set<size_t> facetsL_; // H facet with local G facet

  Array<size_t> vertex_map_;
  Array<size_t> cell_map_;

};

//-----------------------------------------------------------------------------
inline auto MappedManifold::facet_index(Cell const& boundary_cell) -> size_t
{
  dolfin_assert(&boundary_cell.mesh() == this);
  return cell_map_[boundary_cell.index()];
}

//-----------------------------------------------------------------------------
inline auto MappedManifold::vertex_index(Vertex const& boundary_vertex) -> size_t
{
  dolfin_assert(&boundary_vertex.mesh() == this);
  return vertex_map_[boundary_vertex.index()];
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MAPPED_MANIFOLD_H */
