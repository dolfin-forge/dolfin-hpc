// Copyright (C) 2007-2008 Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_DOF_MAP_H
#define __DOLFIN_DOF_MAP_H

#include <dolfin/common/types.h>
#include <dolfin/fem/DofMapCache.h>
#include <dolfin/fem/DofNumbering.h>
#include <dolfin/fem/PeriodicDofsMapping.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/mesh/MeshDependent.h>
#include <dolfin/ufc/ufc.h>

namespace dolfin
{

class Cell;
class Form;
class Mesh;
class SubSytem;
class FiniteElementSpace;

/// This class handles the mapping of degrees of freedom.
/// It wraps a ufc::dofmap on a specific mesh and provides
/// optional precomputation and reordering of dofs.

class DofMap : public ufc::dofmap, public MeshDependent
{

public:
  /// Create dof map on mesh for i-th coefficient of given form
  DofMap( Mesh & mesh, ufc::form const & form, size_t const i );

  /// Create dof map on mesh from UFC object
  /// Ownership of the UFC object is transfered to the instance if the boolean
  /// is set to true, otherwise a clone of the dofmap is created.
  /// In any case the instance the member attribute will be destroyed.
  DofMap( Mesh & mesh, ufc::dofmap & dofmap, bool const owner );

  /// Create dof map on a subspace
  DofMap( DofMap const & dofmap, size_t i );

  /// Create dof map on a subspace for given subsystem
  DofMap( DofMap const &                dofmap,
          std::vector< size_t > const & sub_system,
          size_t &                      offset );

  /// Destructor
  ~DofMap() override;

  /// Check if the element definitions are identical
  auto operator==( DofMap const & other ) const -> bool;
  auto operator!=( DofMap const & other ) const -> bool;

  //--- Instantiation using the dofmap cache

  /// Acquire dofmap from cache for i-th function of the form.
  static auto acquire( Mesh & mesh, Form const & form, size_t const i )
    -> DofMap &;

  /// Acquire dofmap from cache for the given UFC dofmap.
  static auto acquire( Mesh & mesh, ufc::dofmap & dofmap, bool owner )
    -> DofMap &;

  /// Release a token for the given dofmap
  static void release( DofMap & dofmap );

  //--- UFC INTERFACE ---------------------------------------------------------

  /// Return a string identifying the dof map
  auto signature() const -> char const * override;

  /// Return true iff mesh entities of topological dimension d are needed
  auto needs_mesh_entities( size_t d ) const -> bool override;

  /// Return the topological dimension of the associated cell shape
  auto topological_dimension() const -> size_t override;

  /// Return the dimension of the global finite element function space
  auto global_dimension(
    std::vector< size_t > const & num_global_mesh_entities ) const
    -> size_t override;

  /// Return the dimension of the local finite element function space
  /// Return the number of dofs with global support (i.e. global constants)
  auto num_global_support_dofs() const -> size_t override;

  /// Return the dimension of the local finite element function space
  /// for a cell (not including global support dofs)
  auto num_element_support_dofs() const -> size_t override;

  /// Return the dimension of the local finite element function space
  /// for a cell (old version including global support dofs)
  /// was local_dimension() before
  auto num_element_dofs() const -> size_t override;

  /// Return number of facet dofs
  auto num_facet_dofs() const -> size_t override;

  /// Return the number of dofs associated with each cell entity of dimension d
  auto num_entity_dofs( size_t d ) const -> size_t override;

  /// Return the number of dofs associated with the closure
  /// of each cell entity dimension d
  auto num_entity_closure_dofs( size_t d ) const -> size_t override;

  /// Tabulate the local-to-global mapping of dofs on a cell
  void tabulate_dofs( size_t *                      dofs,
                      std::vector< size_t > const & num_global_entities,
                      std::vector< std::vector< size_t > > const &
                        entity_indices ) const override;

  /// Tabulate local-local facet dofs
  void tabulate_facet_dofs( size_t * dofs, size_t local_facet ) const override;

  /// Tabulate the local-to-local mapping of dofs on entity (d, i)
  void tabulate_entity_dofs( size_t * dofs, size_t d, size_t i ) const override;

  // /// Tabulate the coordinates of all dofs on a cell
  // void tabulate_coordinates(real** coordinates,
  //                           const ufc::cell& ufc_cell) const override;

  /// Tabulate the local-to-local mapping of dofs on the closure of entity (d,
  /// i)
  void tabulate_entity_closure_dofs( size_t * dofs,
                                     size_t   d,
                                     size_t   i ) const override;

  //// Return the number of sub dof maps (for a mixed element)
  auto num_sub_dofmaps() const -> size_t override;

  /// Create a new dofmap for sub dof map i (for a mixed element)
  auto create_sub_dofmap( size_t i ) const -> ufc::dofmap * override;

  /// Create a new instance
  auto create() const -> ufc::dofmap * override;

  //--- EXTENSION OF UFC INTERFACE --------------------------------------------

  /// FIXME
  // /// Return the dimension of the local finite element function space
  // size_t macro_local_dimension() const;

  /// Tabulate the local-to-global mapping of dofs on a cell
  /// FIXME
  void tabulate_dofs( size_t *          dofs,
                      ufc::cell const & ufc_cell,
                      Cell const &      cell ) const;

  /// Tabulate the local-to-global mapping of dofs on a cell
  /// FIXME
  void tabulate_dofs( size_t * dofs, UFCCell const & ufc_cell ) const;

  /// Extract sub dof map
  auto create_sub_dofmap( std::vector< size_t > const & sub_system ) const
    -> ufc::dofmap *;

  /// Extract sub dof map and compute the offset local to the reference element
  auto create_sub_dofmap( std::vector< size_t > const & sub_system,
                          size_t & local_offset ) const -> ufc::dofmap *;

  static auto create_sub_dofmap( ufc::dofmap const &           dofmap,
                                 std::vector< size_t > const & sub_system,
                                 size_t & local_offset ) -> ufc::dofmap *;

  /// Get sub dof maps offset (for a mixed element)
  auto sub_dofmaps_dimensions() const -> std::vector< size_t > const &;

  /// Get sub dof maps offset (for a mixed element)
  auto sub_dofmaps_offsets() const -> std::vector< size_t > const &;

  /// Get list of scalar dofmaps ordered by entries
  auto flatten() const -> std::vector< ufc::dofmap const * > const &;

  /// Create flatten representation of given dofmap (append sub dofmaps)
  static void flatten( ufc::dofmap const *                  dofmap,
                       std::vector< ufc::dofmap const * > & stack,
                       size_t                               maxlevel );

  /// Create flatten representation of given dofmap (append sub dofmaps)
  static void flatten( ufc::dofmap const *                  dofmap,
                       std::vector< ufc::dofmap const * > & stack );

  /// Return if the dofmap can be seen as a vector element dofmap
  auto is_vectorizable() const -> bool;

  /// Return if the list of dofmap can be seen as a vector element
  static auto can_vectorize( std::vector< ufc::dofmap const * > flattened )
    -> bool;

  /// Unique identifier
  auto hash() const -> std::string const &;

  /// Return the dofmap local size i.e the process range (only owned dofs)
  auto local_size() const -> size_t;

  //--- Management of local-to-global mapping

  /// Return local to global mapping
  auto dofsmapping() const -> size_t const *;

  /// Return the size of the local to global mapping
  auto dofsmapping_size() const -> size_t;

  //--- Management of periodic dofs mapping

  /// Return local to global mapping
  auto periodic_mapping( FiniteElementSpace const & space ) const
    -> PeriodicDofsMapping const &;

  //--- Dof ownership

  /// Return is the dof is shared, return false if the index is not known (!)
  auto is_shared( size_t index ) const -> bool;

  /// Return is the dof is ghosted, return false if the index is not known (!)
  auto is_ghost( size_t index ) const -> bool;

  //--- Debugging

  /// Return renumbering (used for testing)
  auto get_map() const -> _ordered_map< size_t, size_t >;

  /// Display mapping
  void disp() const;

  /// Return if the dof map has been renumbered
  auto renumbered() const -> bool;

  //---

  /// Returns the dofmap signature corresponding to a given finite element
  static auto make_signature( std::string const & finite_element )
    -> std::string const;

  /// Create unique string identifiers for dofmap from UFC dofmap
  static auto make_hash( Mesh & mesh, ufc::dofmap const & ufc_dofmap )
    -> std::string const;

  //--- Debugging

  /// Check consistency of ghosted entities
  auto check( bool throw_error = false ) -> bool;

private:
  /// Initialize
  void init();

  /// Build dof numbering
  void build();

  // Forward declaration of offset to be update at creation of ufc::dofmap
  size_t offset_;

  // UFC dof map
  ufc::dofmap * const ufc_dofmap_;

  // Use type to allow optimization of dof map ordering
  DofNumbering * const numbering_;

  // Dofmap hash
  std::string const hash_;

  // Sub dof maps offsets
  std::vector< size_t > sub_dofmaps_dims_;

  // Sub dof maps offsets
  std::vector< size_t > sub_dofmaps_offs_;

  //
  mutable std::vector< ufc::dofmap const * > flattened_;

  //
  size_t num_leaf_spaces_;

  // Periodic dofs mapping
  mutable PeriodicDofsMapping * periodic_dofmap_;

  // Provide easy access to map for testing
  _ordered_map< size_t, size_t > map_;
};

//--- INLINES -----------------------------------------------------------------

//-----------------------------------------------------------------------------

inline auto DofMap::operator==( DofMap const & other ) const -> bool
{
  return ( this->hash() == other.hash() );
}

//-----------------------------------------------------------------------------

inline auto DofMap::operator!=( DofMap const & other ) const -> bool
{
  return !( *this == other );
}

//-----------------------------------------------------------------------------

inline auto DofMap::acquire( Mesh & mesh, Form const & form, size_t const i )
  -> DofMap &
{
  return DofMapCache::instance().acquire( mesh, form, i );
}

//-----------------------------------------------------------------------------

inline auto DofMap::acquire( Mesh & mesh, ufc::dofmap & dofmap, bool owner )
  -> DofMap &
{
  return DofMapCache::instance().acquire( mesh, dofmap, owner );
}

//-----------------------------------------------------------------------------

inline void DofMap::release( DofMap & dofmap )
{
  return DofMapCache::instance().release( dofmap );
}

//-----------------------------------------------------------------------------

inline auto DofMap::signature() const -> char const *
{
  return ufc_dofmap_->signature();
}

//-----------------------------------------------------------------------------

inline auto DofMap::needs_mesh_entities( size_t d ) const -> bool
{
  return ufc_dofmap_->needs_mesh_entities( d );
}

//-----------------------------------------------------------------------------

inline auto DofMap::topological_dimension() const -> size_t
{
  return ufc_dofmap_->topological_dimension();
}

//-----------------------------------------------------------------------------

inline auto DofMap::global_dimension(
  std::vector< size_t > const & num_global_mesh_entities ) const -> size_t
{
  return ufc_dofmap_->global_dimension( num_global_mesh_entities );
}

//-----------------------------------------------------------------------------

inline auto DofMap::num_global_support_dofs() const -> size_t
{
  return ufc_dofmap_->num_global_support_dofs();
}

//-----------------------------------------------------------------------------

inline auto DofMap::num_element_support_dofs() const -> size_t
{
  return ufc_dofmap_->num_element_support_dofs();
}

//-----------------------------------------------------------------------------

inline auto DofMap::num_element_dofs() const -> size_t
{
  return ufc_dofmap_->num_element_dofs();
}

//-----------------------------------------------------------------------------

inline auto DofMap::num_facet_dofs() const -> size_t
{
  return ufc_dofmap_->num_facet_dofs();
}

//-----------------------------------------------------------------------------

inline auto DofMap::num_entity_dofs( size_t d ) const -> size_t
{
  return ufc_dofmap_->num_entity_dofs( d );
}

//-----------------------------------------------------------------------------

inline auto DofMap::num_entity_closure_dofs( size_t d ) const -> size_t
{
  return ufc_dofmap_->num_entity_closure_dofs( d );
}

//-----------------------------------------------------------------------------

inline void DofMap::tabulate_dofs(
  size_t *                                     dofs,
  std::vector< size_t > const &                num_global_entities,
  std::vector< std::vector< size_t > > const & entity_indices ) const
{
  ufc_dofmap_->tabulate_dofs( dofs, num_global_entities, entity_indices );
}

//-----------------------------------------------------------------------------

inline void DofMap::tabulate_facet_dofs( size_t * dofs,
                                         size_t   local_facet ) const
{
  ufc_dofmap_->tabulate_facet_dofs( dofs, local_facet );
}

//-----------------------------------------------------------------------------

inline void
  DofMap::tabulate_entity_dofs( size_t * dofs, size_t d, size_t i ) const
{
  ufc_dofmap_->tabulate_entity_dofs( dofs, d, i );
}

//-----------------------------------------------------------------------------

// inline void DofMap::tabulate_coordinates( real **           coordinates,
//                                           const ufc::cell & ufc_cell ) const
// {
//   ufc_dofmap_->tabulate_coordinates( coordinates, ufc_cell );
// }

inline void DofMap::tabulate_entity_closure_dofs( size_t * dofs,
                                                  size_t   d,
                                                  size_t   i ) const
{
  ufc_dofmap_->tabulate_entity_closure_dofs( dofs, d, i );
}

//-----------------------------------------------------------------------------

inline auto DofMap::num_sub_dofmaps() const -> size_t
{
  return ufc_dofmap_->num_sub_dofmaps();
}

//-----------------------------------------------------------------------------

inline auto DofMap::create_sub_dofmap( size_t i ) const -> ufc::dofmap *
{
  if ( ufc_dofmap_-> num_sub_dofmaps() == 0 )
  {
    return ufc_dofmap_->create();
  }
  else
  {
    return ufc_dofmap_->create_sub_dofmap( i );
  }
}

//-----------------------------------------------------------------------------

inline auto DofMap::create() const -> ufc::dofmap *
{
  return ufc_dofmap_->create();
}

//-----------------------------------------------------------------------------

// inline size_t DofMap::macro_local_dimension() const
// {
//   return 2 * ufc_dofmap_->local_dimension();
// }

//-----------------------------------------------------------------------------

inline void DofMap::tabulate_dofs( size_t *          dofs,
                                   ufc::cell const & ufc_cell,
                                   Cell const &      cell ) const
{
  dolfin_assert( dofs != nullptr );
  numbering_->tabulate_dofs( dofs, ufc_cell, cell );
}

//-----------------------------------------------------------------------------

inline void DofMap::tabulate_dofs( size_t *        dofs,
                                   UFCCell const & ufc_cell ) const
{
  dolfin_assert( dofs != nullptr );
  numbering_->tabulate_dofs( dofs, ufc_cell, *ufc_cell );
}

//-----------------------------------------------------------------------------

inline auto DofMap::is_vectorizable() const -> bool
{
  return DofMap::can_vectorize( this->flatten() );
}

//-----------------------------------------------------------------------------

inline auto DofMap::sub_dofmaps_dimensions() const
  -> std::vector< size_t > const &
{
  return sub_dofmaps_dims_;
}

//-----------------------------------------------------------------------------

inline auto DofMap::sub_dofmaps_offsets() const -> std::vector< size_t > const &
{
  return sub_dofmaps_offs_;
}
//-----------------------------------------------------------------------------

inline auto DofMap::local_size() const -> size_t
{
  return numbering_->size();
}

//-----------------------------------------------------------------------------

inline auto DofMap::dofsmapping() const -> size_t const *
{
  return numbering_->block();
}

//-----------------------------------------------------------------------------

inline auto DofMap::dofsmapping_size() const -> size_t
{
  return numbering_->block_size();
}

//-----------------------------------------------------------------------------

inline auto DofMap::get_map() const -> _ordered_map< size_t, size_t >
{
  return map_;
}

//-----------------------------------------------------------------------------

inline auto DofMap::is_shared( size_t index ) const -> bool
{
  return numbering_->is_shared( index );
}

//-----------------------------------------------------------------------------

inline auto DofMap::is_ghost( size_t index ) const -> bool
{
  return numbering_->is_ghost( index );
}

//-----------------------------------------------------------------------------

inline auto DofMap::hash() const -> std::string const &
{
  return hash_;
}

//-----------------------------------------------------------------------------

inline auto DofMap::make_signature( std::string const & finite_element )
  -> std::string const
{
  return "FFC dofmap for " + finite_element;
}

//-----------------------------------------------------------------------------

inline auto DofMap::make_hash( Mesh & mesh, ufc::dofmap const & ufc_dofmap )
  -> std::string const
{
  return std::string( ufc_dofmap.signature() ) + mesh.hash();
}

//-----------------------------------------------------------------------------

}

#endif
