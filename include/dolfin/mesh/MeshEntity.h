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
  ~MeshEntity();

  /// Return mesh associated with mesh entity
  Mesh & mesh();

  /// Return mesh associated with mesh entity
  Mesh const & mesh() const;

  //--- Topology --------------------------------------------------------------

  /// Return topological dimension
  uint dim() const;

  /// Return index of mesh entity
  uint index() const;

  /// Return number of incident mesh entities of given topological dimension
  uint num_entities( uint dim ) const; //!< @tod remove this function

  /// Copy global indices of mesh entities to array
  void get_entities( uint dim, uint * indices ) const;

  /// Copy global indices of mesh entities to array
  void get_entities( uint ** indices ) const;

  /// Return array of indices for incident mesh entities of given topological
  /// dimension
  Array< uint > & entities( uint dim );

  /// Return array of indices for incident mesh entities of given topological
  /// dimension
  Array< uint > const & entities( uint dim ) const;

  /// Check if given entity is incident
  bool incident( MeshEntity const & entity ) const;

  /// Compute local index of given incident entity (-1 if not found)
  int index( MeshEntity const & entity ) const;

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
  uint global_index() const;

  /// Copy global indices of mesh entities to array
  void get_global_entities( uint dim, uint * indices ) const;

  /// Copy global indices of mesh entities to array
  void get_global_entities( uint ** indices ) const;

  /// Return if the mesh entity is owned
  bool is_owned() const;

  /// Return if the mesh entity is shared
  bool is_shared() const;

  /// Return if the mesh entity is ghosted
  bool is_ghost() const;

  /// Return the owner of the mesh entity
  uint owner() const;

  /// Returns a pointer to the adjacent set or null is non-shared
  _set< uint > const * adjacents() const;

  /// Return if the mesh entity has all vertices shared
  bool has_all_vertices_shared() const;

  /// Return if the mesh entity is located on the global mesh boundary
  bool on_boundary() const;

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

inline Mesh & MeshEntity::mesh()
{
  return mesh_;
}

//-----------------------------------------------------------------------------
inline Mesh const & MeshEntity::mesh() const
{
  return mesh_;
}

//-----------------------------------------------------------------------------
inline uint MeshEntity::dim() const
{
  return tdim_;
}

//-----------------------------------------------------------------------------
inline uint MeshEntity::index() const
{
  return index_;
}

//-----------------------------------------------------------------------------
inline uint MeshEntity::num_entities( uint dim ) const
{
  // dolfin_assert(topology_(tdim_, dim).is_initialized());
  dolfin_assert( topology_( tdim_, dim ).min_degree()
                 <= topology_( tdim_, dim ).max_degree() );
  // NOTE: New MeshTopology class auto-creates connectivity on demand.
  return topology_( tdim_, dim ).degree( index_ );
}

//-----------------------------------------------------------------------------
inline Array< uint > & MeshEntity::entities( uint dim )
{
  // dolfin_assert(topology_(tdim_, dim).is_initialized());
  dolfin_assert( topology_( tdim_, dim ).min_degree()
                 <= topology_( tdim_, dim ).max_degree() );
  // NOTE: New MeshTopology class auto-creates connectivity on demand.
  return topology_( tdim_, dim )[index_];
}

//-----------------------------------------------------------------------------
inline Array< uint > const & MeshEntity::entities( uint dim ) const
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
inline bool MeshEntity::incident( MeshEntity const & entity ) const
{
  // Must be in the same mesh to be incident
  if ( &topology_ != &entity.topology_ )
    return false;

  return topology_( tdim_, entity.tdim_ ).incident( index_, entity.index_ );
}

//-----------------------------------------------------------------------------
inline int MeshEntity::index( MeshEntity const & entity ) const
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
inline uint MeshEntity::global_index() const
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
inline bool MeshEntity::is_owned() const
{
  return ( topology_.distributed() ? distdata_[tdim_].is_owned( index_ )
                                   : true );
}

//-----------------------------------------------------------------------------
inline bool MeshEntity::is_shared() const
{
  return ( topology_.distributed() ? distdata_[tdim_].is_shared( index_ )
                                   : false );
}

//-----------------------------------------------------------------------------
inline bool MeshEntity::is_ghost() const
{
  return ( topology_.distributed() ? distdata_[tdim_].is_ghost( index_ )
                                   : false );
}

//-----------------------------------------------------------------------------
inline uint MeshEntity::owner() const
{
  return ( topology_.distributed() ? distdata_[tdim_].get_owner( index_ )
                                   : MPI::rank() );
}

//-----------------------------------------------------------------------------
inline _set< uint > const * MeshEntity::adjacents() const
{
  return ( topology_.distributed() ? distdata_[tdim_].ptr_shared_adj( index_ )
                                   : nullptr );
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_ENTITY_H */
