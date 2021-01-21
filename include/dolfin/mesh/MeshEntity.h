// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MESH_ENTITY_H
#define __DOLFIN_MESH_ENTITY_H

#include <dolfin/common/types.h>
#include <dolfin/mesh/MeshDistributedData.h>
#include <dolfin/mesh/MeshTopology.h>

namespace dolfin
{

class Mesh;

/**
 *
 *  @class  MeshEntity
 *
 *  @brief  A MeshEntity represents a mesh entity associated with a specific
 *          topological dimension of some mesh.
 *
 */

class MeshEntity
{

public:
  /// Constructor
  MeshEntity( Mesh & mesh, uint dim, uint index );

  /// Destructor
  ~MeshEntity() = default;

  /// Return mesh associated with mesh entity
  auto mesh() -> Mesh &;

  /// Return mesh associated with mesh entity
  auto mesh() const -> Mesh const &;

  //--- Topology --------------------------------------------------------------

  /// Return topological dimension
  auto dim() const -> uint;

  /// Return index of mesh entity
  auto index() const -> uint;

  /// Return number of incident mesh entities of given topological dimension
  auto num_entities( uint dim ) const -> uint; //!< @tod remove this function

  /// Copy global indices of mesh entities to array
  void get_entities( uint dim, uint * indices ) const;

  /// Copy global indices of mesh entities to array
  void get_entities( uint ** indices ) const;

  /// Return array of indices for incident mesh entities of given topological
  /// dimension
  auto entities( uint dim ) -> Array< uint > &;

  /// Return array of indices for incident mesh entities of given topological
  /// dimension
  auto entities( uint dim ) const -> Array< uint > const &;

  /// Check if given entity is incident
  auto incident( MeshEntity const & entity ) const -> bool;

  /// Compute local index of given incident entity (-1 if not found)
  auto index( MeshEntity const & entity ) const -> int;

  //--- Geometry --------------------------------------------------------------

  /*
   * TBD: add convenience function for accessing vertex coordinates
   *
   */

  //---------------------------------------------------------------------------

  /*
   * Convenience functions which can be called in serial and parallel.
   * For code portions within which the mesh is known to be distributed, calling
   * functions from the mesh distributed data is recommended.
   *
   */

  /// Return global index of mesh entity
  auto global_index() const -> uint;

  /// Copy global indices of mesh entities to array
  void get_global_entities( uint dim, uint * indices ) const;

  /// Copy global indices of mesh entities to array
  void get_global_entities( uint ** indices ) const;

  /// Return if the mesh entity is owned
  auto is_owned() const -> bool;

  /// Return if the mesh entity is shared
  auto is_shared() const -> bool;

  /// Return if the mesh entity is ghosted
  auto is_ghost() const -> bool;

  /// Return the owner of the mesh entity
  auto owner() const -> uint;

  /// Returns a pointer to the adjacent set or null is non-shared
  auto adjacents() const -> _set< uint > const *;

  /// Return if the mesh entity has all vertices shared
  auto has_all_vertices_shared() const -> bool;

  /// Return if the mesh entity is located on the global mesh boundary
  auto on_boundary() const -> bool;

  //---------------------------------------------------------------------------

  ///
  void disp() const;

protected:
  // Friends
  friend class MeshEntityIterator;

  // Mesh
  Mesh & mesh_;

  // Mesh topology
  MeshTopology & topology_;

  // Topological dimension
  uint const tdim_;

  // Geometric dimension
  uint const gdim_;

  // Pointer to mesh distributed data if applicable
  MeshDistributedData const & distdata_;

  // Index of entity within topological dimension
  uint index_;
};

//--- INLINES -----------------------------------------------------------------

inline auto MeshEntity::mesh() -> Mesh &
{
  return mesh_;
}

//-----------------------------------------------------------------------------
inline auto MeshEntity::mesh() const -> Mesh const &
{
  return mesh_;
}

//-----------------------------------------------------------------------------
inline auto MeshEntity::dim() const -> uint
{
  return tdim_;
}

//-----------------------------------------------------------------------------
inline auto MeshEntity::index() const -> uint
{
  return index_;
}

//-----------------------------------------------------------------------------
inline auto MeshEntity::num_entities( uint dim ) const -> uint
{
  // dolfin_assert(topology_(tdim_, dim).is_initialized());
  dolfin_assert( topology_( tdim_, dim ).min_degree()
                 <= topology_( tdim_, dim ).max_degree() );
  // NOTE: New MeshTopology class auto-creates connectivity on demand.
  return topology_( tdim_, dim ).degree( index_ );
}

//-----------------------------------------------------------------------------
inline auto MeshEntity::entities( uint dim ) -> Array< uint > &
{
  // dolfin_assert(topology_(tdim_, dim).is_initialized());
  dolfin_assert( topology_( tdim_, dim ).min_degree()
                 <= topology_( tdim_, dim ).max_degree() );
  // NOTE: New MeshTopology class auto-creates connectivity on demand.
  return topology_( tdim_, dim )[index_];
}

//-----------------------------------------------------------------------------
inline auto MeshEntity::entities( uint dim ) const -> Array< uint > const &
{
  // dolfin_assert(topology_(tdim_, dim).is_initialized());
  dolfin_assert( topology_( tdim_, dim ).min_degree()
                 <= topology_( tdim_, dim ).max_degree() );
  // NOTE: New MeshTopology class auto-creates connectivity on demand.
  return topology_( tdim_, dim )[index_];
}

//-----------------------------------------------------------------------------
inline void MeshEntity::get_entities( uint dim, uint * indices ) const
{
  Array< uint > const & e = topology_( tdim_, dim )[index_];
  dolfin_assert( indices != nullptr );
  std::copy( e.begin(), e.end(), indices );
}

//-----------------------------------------------------------------------------
inline void MeshEntity::get_entities( uint ** indices ) const
{
  dolfin_assert( indices != nullptr );
  for ( uint dim = 0; dim < tdim_; ++dim )
  {
    Array< Array< uint > > const & e = topology_( tdim_, dim )();
    dolfin_assert( indices[dim] != nullptr );
    std::copy( e[index_].begin(), e[index_].end(), indices[dim] );
  }
  indices[tdim_][0] = index_;
}

//-----------------------------------------------------------------------------
inline auto MeshEntity::incident( MeshEntity const & entity ) const -> bool
{
  // Must be in the same mesh to be incident
  if ( &topology_ != &entity.topology_ )
    return false;

  return topology_( tdim_, entity.tdim_ ).incident( index_, entity.index_ );
}

//-----------------------------------------------------------------------------
inline auto MeshEntity::index( MeshEntity const & entity ) const -> int
{
  // Must be in the same mesh to be incident
  if ( &topology_ != &entity.topology_ )
  {
    error(
      "Unable to compute index of an entity defined on a different mesh." );
  }

  return topology_( tdim_, entity.tdim_ ).index( index_, entity.index_ );
}

//-----------------------------------------------------------------------------
inline auto MeshEntity::global_index() const -> uint
{
  return ( topology_.distributed() ? distdata_[tdim_].get_global( index_ )
                                   : index_ );
}

//-----------------------------------------------------------------------------
inline void MeshEntity::get_global_entities( uint dim, uint * indices ) const
{
  // Get list of entities for given topological dimension
  if ( topology_.distributed() )
  {
    Connectivity const & mc = topology_( tdim_, dim );
    distdata_[dim].get_global(
      mc.degree( index_ ), mc[index_].data(), indices );
  }
  else
  {
    get_entities( dim, indices );
  }
}

//-----------------------------------------------------------------------------
inline void MeshEntity::get_global_entities( uint ** indices ) const
{
  // Get list of entities for given topological dimension
  if ( topology_.distributed() )
  {
    for ( uint d = 0; d < tdim_; ++d )
    {
      Connectivity const & mc = topology_( tdim_, d );
      distdata_[d].get_global(
        mc.degree( index_ ), mc[index_].data(), indices[d] );
    }
    indices[tdim_][0] = distdata_[tdim_].get_global( index_ );
  }
  else
  {
    get_entities( indices );
  }
}

//-----------------------------------------------------------------------------
inline auto MeshEntity::is_owned() const -> bool
{
  return ( topology_.distributed() ? distdata_[tdim_].is_owned( index_ )
                                   : true );
}

//-----------------------------------------------------------------------------
inline auto MeshEntity::is_shared() const -> bool
{
  return ( topology_.distributed() ? distdata_[tdim_].is_shared( index_ )
                                   : false );
}

//-----------------------------------------------------------------------------
inline auto MeshEntity::is_ghost() const -> bool
{
  return ( topology_.distributed() ? distdata_[tdim_].is_ghost( index_ )
                                   : false );
}

//-----------------------------------------------------------------------------
inline auto MeshEntity::owner() const -> uint
{
  return ( topology_.distributed() ? distdata_[tdim_].get_owner( index_ )
                                   : MPI::rank() );
}

//-----------------------------------------------------------------------------
inline auto MeshEntity::adjacents() const -> _set< uint > const *
{
  return ( topology_.distributed() ? distdata_[tdim_].ptr_shared_adj( index_ )
                                   : nullptr );
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_ENTITY_H */
