// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

// This class has been extensively tortured since 2015 due to flawed design and
// non-robust implementation.
// Additionally the topology computation was rewritten to respect abstraction
// layers and avoid storing cell - cell neighbours connectivities.
// The updated interface gives access to the distributed data to allow writing
// generic serial/parallel code with no requirement of parallel guards for
// sections reserved to distributed meshes.
// Additionally it benefits from several inconsistency fixes merged from mesh
// distributed data.

#ifndef __DOLFIN_MESH_TOPOLOGY_H
#define __DOLFIN_MESH_TOPOLOGY_H

#include <dolfin/common/Clonable.h>
#include <dolfin/common/Distributed.h>
#include <dolfin/common/Tokenized.h>
#include <dolfin/common/types.h>
#include <dolfin/log/log.h>
#include <dolfin/mesh/Connectivity.h>
#include <dolfin/mesh/MeshDistributedData.h>
#include <dolfin/mesh/celltypes/CellType.h>

namespace dolfin
{

class Mesh;
class Connectivity;

//-----------------------------------------------------------------------------

/**
 *
 *  @class  MeshTopology
 *
 *  @brief  MeshTopology stores the topology of a mesh, consisting of mesh
 *          entities and connectivity (incidence relations of mesh entities).
 *          Note that the mesh entities do not need to be stored, only the
 *          number of entities and the connectivity. Any numbering scheme for
 *          the mesh entities is stored separately in a MeshFunction over the
 *          entities.
 *          A mesh entity e may be identified globally as a pair e = (dim, i),
 *          where dim is the topological dimension and i is the index of the
 *          entity within that topological dimension.
 *
 */

class MeshTopology : public Clonable< MeshTopology >,
                     public Distributed< MeshTopology >
{
  // Set max connectivity matrix dimension
  static constexpr size_t CMAX = 4;

  // Save some limbo at MeshEntity construction until classes are rewritten
  friend class MeshEntity;

public:
  /// Create mesh topology for given cell type
  MeshTopology( CellType const & type, Comm & comm, bool frozen );

  /// Copy constructor
  MeshTopology( MeshTopology const & other );

  /// Destructor
  ~MeshTopology() override;

  /// Swap instances
  friend auto swap( MeshTopology & a, MeshTopology & b ) -> void;

  /// Assignment (Disabled)
  auto operator=( const MeshTopology & other ) -> MeshTopology &;

  /// Equality
  auto operator==( MeshTopology const & other ) const -> bool;

  /// Non-equality
  auto operator!=( MeshTopology const & other ) const -> bool;

  /// Set topology entities for given dimension
  /// Optionally specify the global number of entities for a distributed mesh.
  /// If the topology is not distributed, any value different than zero or the
  /// number of local entities will trigger an error.
  auto init( size_t dim, size_t nlocal, size_t nglobal = 0 ) -> void;

  /// Finalize: check, reorder, renumber
  auto finalize() -> void;

  ///
  auto type() const -> CellType const &;

  /// Remap local entities of given dimension
  auto remap( size_t d0, std::vector< size_t > const & mapping ) -> void;

  //--- Connectivity ----------------------------------------------------------

  /// Return connectivity for given pair of topological dimensions
  auto operator()( size_t d0, size_t d1 ) -> Connectivity &;

  /// Return connectivity for given pair of topological dimensions
  auto operator()( size_t d0, size_t d1 ) const -> Connectivity const &;

  /// Return topological dimension
  auto dim() const -> size_t;

  /// number of entities in all topological dimensions
  auto num_entities() const -> std::vector< size_t > const &;

  /// update number of entities in all topological dimensions
  /// FIXME this shouldnt be const, but this class is somehow built to not care
  /// about that...
  auto num_entities_update() const -> void;

  /// Return number of entities in the local topology for given dimension
  auto size( size_t dim ) const -> size_t;

  /// Return pointer to connectivity for given pair
  auto connectivity( size_t d0, size_t d1 = 0 ) -> Connectivity *;

  /// Return pointer to connectivity for given pair (const)
  auto connectivity( size_t d0, size_t d1 = 0 ) const
    -> Connectivity const *;

  //--- Distributed data ------------------------------------------------------

  /// Return mesh distribution data if the topology is distributed
  auto distdata() -> MeshDistributedData &;

  /// Return mesh distribution data if the topology is distributed (const)
  auto distdata() const -> MeshDistributedData const &;

  /// Return number of entities in the global topology for given dimension
  auto global_size( size_t dim ) const -> size_t;

  /// Return offset of global indices on current rank
  auto offset( size_t dim ) const -> size_t;

  /// Return number of given entities
  auto num_owned( size_t dim ) const -> size_t;

  /// Return number of given entities
  auto num_shared( size_t dim ) const -> size_t;

  /// Return number of given entities
  auto num_ghost( size_t dim ) const -> size_t;

  //---------------------------------------------------------------------------

  /// Display data
  auto disp() const -> void;

  //--- TOKENIZED -------------------------------------------------------------

  /// Return token identifying the internal state of mesh topology
  auto token() const -> int;

private:
  /// Update token value
  auto update_token() -> void;

  //---------------------------------------------------------------------------

  /// Compute connectivity for given pair of topological dimensions
  auto compute( size_t d0, size_t d1 ) const -> Connectivity const *;

  /// Compute entities for given topological dimension
  auto entities( size_t di ) const -> Connectivity const *;

  /// Compute transpose for given pair of topological dimensions
  auto transpose( size_t d0, size_t d1 ) const -> Connectivity const *;

  /// Compute connectivity for given triple of topological dimensions
  auto intersection( size_t d0, size_t di, size_t d1 ) const
    -> Connectivity const *;

public:
  /// Force renumbering of mesh topology entities
  /// @todo public for the moment but this just legacy of bad design
  /// FIXME this shouldnt be const, but this class is somehow built to not care
  /// about that...
  auto renumber() const -> void;

private:
  /// Force reordering of mesh topology connectivities
  /// FIXME this shouldnt be const, but this class is somehow built to not care
  /// about that...
  auto reorder() const -> void;

  ///
  CellType const * type_;

  /// Topological dimension
  size_t dim_;

  // Topology cannot be modified
  bool frozen_;

  /// number of entities in each (connectivity) dimension
  /// FIXME this shouldnt be mutable, but this class is somehow built to not care
  /// about that...
  mutable std::vector< size_t > num_entities_;

  /// Connectivity for pairs of topological dimensions
  mutable Connectivity * C_[CMAX][CMAX];

  /// Distributed mesh topology data
  MeshDistributedData distdata_;

  //
  int timestamp_;
};

//-----------------------------------------------------------------------------

inline auto MeshTopology::type() const -> CellType const &
{
  dolfin_assert( type_ );
  return *type_;
}

//-----------------------------------------------------------------------------

inline auto MeshTopology::operator()( size_t d0, size_t d1 ) -> Connectivity &
{
  dolfin_assert( d0 <= dim_ && d1 <= dim_ );
  if ( !connectivity( d0, d1 ) )
  {
    compute( d0, d1 );
  }
  return *connectivity( d0, d1 );
}

//-----------------------------------------------------------------------------

inline auto MeshTopology::operator()( size_t d0, size_t d1 ) const
  -> Connectivity const &
{
  dolfin_assert( d0 <= dim_ && d1 <= dim_ );
  return *compute( d0, d1 );
}

//-----------------------------------------------------------------------------

inline auto MeshTopology::dim() const -> size_t
{
  return dim_;
}

//-----------------------------------------------------------------------------

inline auto MeshTopology::num_entities() const -> std::vector< size_t > const &
{
  dolfin_assert( num_entities_.size() == dim_ + 1 );
  return num_entities_;
}

//-----------------------------------------------------------------------------

inline auto MeshTopology::size( size_t dim ) const -> size_t
{
  dolfin_assert( dim <= dim_ );
  return ( *this )( dim, 0 ).order();
}

//-----------------------------------------------------------------------------

inline auto MeshTopology::connectivity( size_t d0, size_t d1 ) -> Connectivity *
{
  dolfin_assert( d0 <= dim_ );
  dolfin_assert( d1 <= dim_ );
  return C_[d0][d1];
}

//-----------------------------------------------------------------------------

inline auto MeshTopology::connectivity( size_t d0, size_t d1 ) const
  -> Connectivity const *
{
  dolfin_assert( d0 <= dim_ );
  dolfin_assert( d1 <= dim_ );
  return C_[d0][d1];
}

//-----------------------------------------------------------------------------

inline auto MeshTopology::distdata() -> MeshDistributedData &
{
  if ( not distributed() )
  {
    error( "MeshTopology : returning distributed data of serial mesh" );
  }
  return distdata_;
}

//-----------------------------------------------------------------------------

inline auto MeshTopology::distdata() const -> MeshDistributedData const &
{
  if ( not distributed() )
  {
    error( "MeshTopology : returning distributed data of serial mesh" );
  }
  return distdata_;
}

//-----------------------------------------------------------------------------

inline auto MeshTopology::global_size( size_t dim ) const -> size_t
{
  return ( distributed() ? distdata_[dim].global_size() : this->size( dim ) );
}

//-----------------------------------------------------------------------------

inline auto MeshTopology::offset( size_t dim ) const -> size_t
{
  return ( distributed() ? distdata_[dim].offset() : 0 );
}

//-----------------------------------------------------------------------------

inline auto MeshTopology::num_owned( size_t dim ) const -> size_t
{
  return ( distributed() ? distdata_[dim].num_owned() : this->size( dim ) );
}

//-----------------------------------------------------------------------------

inline auto MeshTopology::num_shared( size_t dim ) const -> size_t
{
  return ( distributed() ? distdata_[dim].num_shared() : 0 );
}

//-----------------------------------------------------------------------------

inline auto MeshTopology::num_ghost( size_t dim ) const -> size_t
{
  return ( distributed() ? distdata_[dim].num_ghost() : 0 );
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_TOPOLOGY_H */
